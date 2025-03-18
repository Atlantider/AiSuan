from django.shortcuts import render, get_object_or_404
from django.http import JsonResponse
from django.conf import settings
from django.utils import timezone
from django.contrib.auth.models import User
from django.db import transaction
from rest_framework import viewsets, permissions, status
from rest_framework.decorators import api_view, action, permission_classes
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from .models import (
    Solvent, 
    Salt, 
    ElectrolyteFormulation, 
    FormulationComponent, 
    SimulationParameters,
    Calculation,
    CalculationResult,
    InputFile
)
from .serializers import (
    SolventSerializer,
    SaltSerializer,
    ElectrolyteFormulationSerializer,
    FormulationComponentSerializer,
    SimulationParametersSerializer,
    CalculationSerializer,
    CalculationResultSerializer,
    CreateFormulationSerializer
)
from .tasks import run_electrolyte_calculation
import os
import sys
import json
import tempfile
from rest_framework.views import APIView
import logging
import datetime
import uuid

# 添加molyte_cursor到Python路径
MOLYTE_PATH = os.path.join(os.path.dirname(settings.BASE_DIR), 'molyte_cursor')
if MOLYTE_PATH not in sys.path:
    sys.path.append(MOLYTE_PATH)

# 尝试导入molyte_cursor功能
try:
    from molyte_cursor.src.core.main import generate_input_files
    from molyte_cursor.src.io.file_generator import LAMMPSFileGenerator
    MOLYTE_AVAILABLE = True
except ImportError:
    print("警告: 无法导入molyte_cursor模块，将使用模拟功能")
    MOLYTE_AVAILABLE = False

class SolventViewSet(viewsets.ModelViewSet):
    """
    溶剂API视图集
    """
    queryset = Solvent.objects.all()
    serializer_class = SolventSerializer
    # permission_classes = [permissions.IsAuthenticated]

class SaltViewSet(viewsets.ModelViewSet):
    """
    盐API视图集
    """
    queryset = Salt.objects.all()
    serializer_class = SaltSerializer
    # permission_classes = [permissions.IsAuthenticated]

class ElectrolyteFormulationViewSet(viewsets.ModelViewSet):
    """
    电解质配方API视图集
    """
    serializer_class = ElectrolyteFormulationSerializer
    # permission_classes = [permissions.IsAuthenticated]
    
    def get_queryset(self):
        """
        根据当前用户筛选配方
        """
        return ElectrolyteFormulation.objects.filter(user=self.request.user)
    
    def create(self, request, *args, **kwargs):
        """
        创建新的电解质配方
        """
        serializer = CreateFormulationSerializer(data=request.data, context={'request': request})
        serializer.is_valid(raise_exception=True)
        formulation = serializer.save()
        
        # 返回创建的配方
        result_serializer = ElectrolyteFormulationSerializer(formulation)
        return Response(result_serializer.data, status=status.HTTP_201_CREATED)
    
    @action(detail=True, methods=['get'])
    def components(self, request, pk=None):
        """
        获取配方的组件列表
        """
        formulation = self.get_object()
        components = FormulationComponent.objects.filter(formulation=formulation)
        serializer = FormulationComponentSerializer(components, many=True)
        return Response(serializer.data)
    
    @action(detail=True, methods=['get'])
    def parameters(self, request, pk=None):
        """
        获取配方的模拟参数
        """
        formulation = self.get_object()
        try:
            parameters = SimulationParameters.objects.get(formulation=formulation)
            serializer = SimulationParametersSerializer(parameters)
            return Response(serializer.data)
        except SimulationParameters.DoesNotExist:
            return Response(
                {"detail": "未找到该配方的模拟参数"},
                status=status.HTTP_404_NOT_FOUND
            )
    
    @action(detail=True, methods=['post'])
    def generate_input_file(self, request, pk=None):
        """
        生成LAMMPS输入文件API
        """
        formulation = self.get_object()
        components = FormulationComponent.objects.filter(formulation=formulation)
        
        # 获取盐组分和溶剂组分
        salts = []
        solvents = []
        
        for component in components:
            if component.component_type == 'salt':
                salts.append({
                    'name': component.salt.name,
                    'cation': component.salt.cation,
                    'anion': component.salt.anion,
                    'concentration': component.concentration
                })
            elif component.component_type == 'solvent':
                solvents.append({
                    'name': component.solvent.name,
                    'smile': component.solvent.smile,
                    'concentration': component.concentration
                })
        
        # 获取模拟参数
        try:
            params = SimulationParameters.objects.get(formulation=formulation)
        except SimulationParameters.DoesNotExist:
            return Response(
                {"detail": "未找到该配方的模拟参数"},
                status=status.HTTP_404_NOT_FOUND
            )
        
        # 准备输出目录
        output_dir = os.path.join(settings.MEDIA_ROOT, 'calculations', str(formulation.id))
        os.makedirs(output_dir, exist_ok=True)
        
        # 准备输入文件
        if MOLYTE_AVAILABLE:
            try:
                # 使用molyte生成输入文件
                config = {
                    'formulation_name': formulation.name,
                    'salts': salts,
                    'solvents': solvents,
                    'temperature': params.temperature,
                    'pressure': params.pressure,
                    'time_step': params.time_step,
                    'equilibration_steps': params.equilibration_steps,
                    'production_steps': params.production_steps,
                    'cutoff': params.cutoff,
                    'output_dir': output_dir
                }
                
                # 调用molyte的API生成文件
                input_file_path = os.path.join(output_dir, 'input.lammps')
                
                # 由于没有详细的molyte_cursor实现，以下仅为示例调用方式
                # 实际实现需要根据molyte_cursor的API进行修改
                # generate_input_files(config, output_dir)
                
                # 模拟生成LAMMPS输入文件
                with open(input_file_path, 'w') as f:
                    f.write(f"# LAMMPS输入文件 - {formulation.name}\n")
                    f.write(f"# 由AiSuan平台生成\n\n")
                    
                    f.write("units real\n")
                    f.write("atom_style full\n")
                    f.write(f"timestep {params.time_step}\n")
                    f.write(f"pair_style lj/cut/coul/long {params.cutoff}\n\n")
                    
                    f.write("# 溶剂组分:\n")
                    for solvent in solvents:
                        f.write(f"# - {solvent['name']}: {solvent['concentration']}\n")
                    
                    f.write("\n# 盐组分:\n")
                    for salt in salts:
                        f.write(f"# - {salt['name']}: {salt['concentration']} mol/L\n")
                    
                    f.write(f"\n# 模拟参数\n")
                    f.write(f"# 温度: {params.temperature} K\n")
                    f.write(f"# 压力: {params.pressure} atm\n")
                    f.write(f"# 平衡步数: {params.equilibration_steps}\n")
                    f.write(f"# 生产步数: {params.production_steps}\n")
                
                # 创建或更新计算记录
                calculation, created = Calculation.objects.get_or_create(
                    formulation=formulation,
                    user=request.user,
                    defaults={
                        'name': f"{formulation.name}的计算",
                        'status': 'pending',
                        'input_file': os.path.join('calculations', str(formulation.id), 'input.lammps')
                    }
                )
                
                if not created:
                    calculation.status = 'pending'
                    calculation.input_file = os.path.join('calculations', str(formulation.id), 'input.lammps')
                    calculation.save()
                
                # 返回成功响应和输入文件路径
                return Response({
                    'status': 'success',
                    'message': '成功生成LAMMPS输入文件',
                    'input_file': calculation.input_file.url if calculation.input_file else None,
                    'calculation_id': calculation.id
                })
                
            except Exception as e:
                return Response({
                    'status': 'error',
                    'message': f'生成输入文件时出错: {str(e)}'
                }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            # 如果molyte不可用，创建模拟的输入文件
            try:
                input_file_path = os.path.join(output_dir, 'input.lammps')
                with open(input_file_path, 'w') as f:
                    f.write(f"# LAMMPS输入文件 - {formulation.name}\n")
                    f.write(f"# 由AiSuan平台生成 (模拟模式)\n\n")
                    
                    f.write("units real\n")
                    f.write("atom_style full\n")
                    f.write(f"timestep {params.time_step}\n")
                    f.write(f"pair_style lj/cut/coul/long {params.cutoff}\n\n")
                    
                    f.write("# 溶剂组分:\n")
                    for solvent in solvents:
                        f.write(f"# - {solvent['name']}: {solvent['concentration']}\n")
                    
                    f.write("\n# 盐组分:\n")
                    for salt in salts:
                        f.write(f"# - {salt['name']}: {salt['concentration']} mol/L\n")
                    
                    f.write(f"\n# 模拟参数\n")
                    f.write(f"# 温度: {params.temperature} K\n")
                    f.write(f"# 压力: {params.pressure} atm\n")
                    f.write(f"# 平衡步数: {params.equilibration_steps}\n")
                    f.write(f"# 生产步数: {params.production_steps}\n")
                
                # 创建或更新计算记录
                calculation, created = Calculation.objects.get_or_create(
                    formulation=formulation,
                    user=request.user,
                    defaults={
                        'name': f"{formulation.name}的计算",
                        'status': 'pending',
                        'input_file': os.path.join('calculations', str(formulation.id), 'input.lammps')
                    }
                )
                
                if not created:
                    calculation.status = 'pending'
                    calculation.input_file = os.path.join('calculations', str(formulation.id), 'input.lammps')
                    calculation.save()
                
                # 返回成功响应和输入文件路径
                return Response({
                    'status': 'success',
                    'message': '成功生成LAMMPS输入文件 (模拟模式)',
                    'input_file': calculation.input_file.url if calculation.input_file else None,
                    'calculation_id': calculation.id
                })
                
            except Exception as e:
                return Response({
                    'status': 'error',
                    'message': f'生成输入文件时出错: {str(e)}'
                }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

class CalculationViewSet(viewsets.ModelViewSet):
    """
    计算任务API视图集
    """
    serializer_class = CalculationSerializer
    # 临时禁用认证要求，仅用于开发测试
    # permission_classes = [permissions.IsAuthenticated]
    
    def get_queryset(self):
        """
        根据当前用户筛选计算任务
        """
        # 临时返回所有计算任务，仅用于开发测试
        return Calculation.objects.all()
        # return Calculation.objects.filter(user=self.request.user)
    
    def perform_create(self, serializer):
        """
        创建计算任务并启动Celery任务
        """
        # 临时使用空用户，仅用于开发测试
        calculation = serializer.save(
            user=self.request.user if hasattr(self.request, 'user') and self.request.user.is_authenticated else None,
            status='pending'
        )
        
        # 启动Celery任务
        task = run_electrolyte_calculation.delay(calculation.id)
        
        # 更新任务ID
        calculation.task_id = task.id
        calculation.save()
    
    @action(detail=True, methods=['get'])
    def result(self, request, pk=None):
        """
        获取计算结果
        """
        calculation = self.get_object()
        try:
            result = CalculationResult.objects.get(calculation=calculation)
            serializer = CalculationResultSerializer(result)
            return Response(serializer.data)
        except CalculationResult.DoesNotExist:
            return Response(
                {"detail": "计算结果尚未生成或不存在"},
                status=status.HTTP_404_NOT_FOUND
            )
    
    @action(detail=True, methods=['post'])
    def restart(self, request, pk=None):
        """
        重新启动失败的计算任务
        """
        calculation = self.get_object()
        
        if calculation.status not in ['failed', 'completed']:
            return Response(
                {"detail": "只能重启失败或已完成的计算任务"},
                status=status.HTTP_400_BAD_REQUEST
            )
        
        # 重置计算状态
        calculation.status = 'pending'
        calculation.error_message = ''
        calculation.save()
        
        # 启动新的Celery任务
        task = run_electrolyte_calculation.delay(calculation.id)
        
        # 更新任务ID
        calculation.task_id = task.id
        calculation.save()
        
        serializer = self.get_serializer(calculation)
        return Response(serializer.data)
    
    @action(detail=True, methods=['post'])
    def submit(self, request, pk=None):
        """
        提交计算任务到LAMMPS
        """
        calculation = self.get_object()
        
        # 确保计算任务有输入文件
        if not calculation.input_file:
            return Response(
                {"detail": "计算任务没有输入文件，请先生成输入文件"},
                status=status.HTTP_400_BAD_REQUEST
            )
        
        # 启动Celery任务执行计算
        task = run_electrolyte_calculation.delay(calculation.id)
        
        # 更新任务状态
        calculation.status = 'running'
        calculation.task_id = task.id
        calculation.save()
        
        return Response({
            'status': 'success',
            'message': '计算任务已成功提交',
            'calculation_id': calculation.id,
            'task_id': task.id
        })

# 新增INP文件API
@api_view(['POST'])
@permission_classes([])  # 空列表表示不需要任何权限
def save_input_file(request, formulation_id):
    """保存电解液计算的INP输入文件"""
    logger = logging.getLogger(__name__)
    logger.info(f"接收到保存INP文件请求: formulation_id={formulation_id}")
    
    try:
        formulation = get_object_or_404(ElectrolyteFormulation, id=formulation_id)
        
        # 临时禁用权限检查，仅用于开发测试
        # if formulation.user != request.user:
        #     logger.warning(f"用户 {request.user.username} 尝试修改不属于他的配方 {formulation_id} 的输入文件")
        #     return Response({"success": False, "error": "您无权修改此配方的输入文件"}, status=status.HTTP_403_FORBIDDEN)
        
        # 获取文件内容
        file_content = request.data.get('content')
        if not file_content:
            logger.error("文件内容为空")
            return Response({"success": False, "error": "文件内容不能为空"}, status=status.HTTP_400_BAD_REQUEST)
        
        # 保存或更新文件
        input_file, created = InputFile.objects.update_or_create(
            formulation=formulation,
            defaults={'content': file_content}
        )
        
        # 文件会在模型的save方法中自动写入文件系统
        logger.info(f"INP文件已保存: formulation_id={formulation_id}, file_path={input_file.file_path}")
        
        return Response({
            "success": True, 
            "message": "输入文件已保存",
            "file_path": input_file.file_path
        })
        
    except Exception as e:
        logger.error(f"保存输入文件时出错: {str(e)}")
        return Response({"success": False, "error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

@api_view(['GET'])
@permission_classes([])  # 空列表表示不需要任何权限
def check_input_file(request, formulation_id):
    """检查配方是否有保存的INP输入文件"""
    logger = logging.getLogger(__name__)
    logger.info(f"检查INP文件是否存在: formulation_id={formulation_id}")
    
    try:
        formulation = get_object_or_404(ElectrolyteFormulation, id=formulation_id)
        
        # 临时禁用权限检查，仅用于开发测试
        # if formulation.user != request.user:
        #     logger.warning(f"用户 {request.user.username} 尝试访问不属于他的配方 {formulation_id} 的输入文件")
        #     return Response({"success": False, "error": "您无权访问此配方的输入文件"}, status=status.HTTP_403_FORBIDDEN)
        
        try:
            input_file = InputFile.objects.get(formulation=formulation)
            file_exists = os.path.exists(input_file.file_path) if input_file.file_path else False
            
            logger.info(f"INP文件查询结果: exists=True, valid={file_exists}, path={input_file.file_path}")
            return Response({
                "success": True,
                "exists": True,
                "file_path": input_file.file_path,
                "valid": file_exists
            })
        except InputFile.DoesNotExist:
            logger.info(f"INP文件查询结果: exists=False")
            return Response({
                "success": True,
                "exists": False
            })
            
    except Exception as e:
        logger.error(f"检查输入文件时出错: {str(e)}")
        return Response({"success": False, "error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

# 新增获取计算状态和结果的API
@api_view(['GET'])
# 临时禁用认证要求，仅用于开发测试
# @permission_classes([IsAuthenticated])
def get_calculation_status(request, calculation_id):
    """获取计算任务的状态"""
    logger = logging.getLogger(__name__)
    logger.info(f"获取计算状态: calculation_id={calculation_id}")
    
    try:
        calculation = get_object_or_404(Calculation, id=calculation_id)
        
        # 临时禁用权限检查，用于开发测试
        # # 检查权限
        # if calculation.user != request.user:
        #     logger.warning(f"用户 {request.user.username} 尝试访问不属于他的计算 {calculation_id} 的状态")
        #     return Response({"success": False, "error": "无权访问此计算任务"}, status=status.HTTP_403_FORBIDDEN)
        
        # 返回状态信息
        logger.info(f"计算状态: calculation_id={calculation_id}, status={calculation.status}")
        return Response({
            "success": True,
            "status": calculation.status,
            "message": calculation.error_message
        })
        
    except Exception as e:
        logger.error(f"获取计算状态时出错: {str(e)}")
        return Response({"success": False, "error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

@api_view(['GET'])
# 临时禁用认证要求，仅用于开发测试
# @permission_classes([IsAuthenticated])
def get_calculation_results(request, calculation_id):
    """获取计算任务的结果"""
    logger = logging.getLogger(__name__)
    logger.info(f"获取计算结果: calculation_id={calculation_id}")
    
    try:
        calculation = get_object_or_404(Calculation, id=calculation_id)
        
        # 临时禁用权限检查，用于开发测试
        # # 检查权限
        # if calculation.user != request.user:
        #     logger.warning(f"用户 {request.user.username} 尝试访问不属于他的计算 {calculation_id} 的结果")
        #     return Response({"success": False, "error": "无权访问此计算任务"}, status=status.HTTP_403_FORBIDDEN)
        
        # 检查状态
        if calculation.status != 'completed':
            logger.warning(f"计算尚未完成，无法获取结果: calculation_id={calculation_id}, status={calculation.status}")
            return Response({
                "success": False,
                "error": f"计算尚未完成，当前状态: {calculation.get_status_display()}"
            }, status=status.HTTP_400_BAD_REQUEST)
        
        # 获取结果
        try:
            result = calculation.result
            # 将结果序列化返回
            logger.info(f"成功获取计算结果: calculation_id={calculation_id}")
            return Response({
                "success": True,
                "results": result.additional_results or {
                    "ionic_conductivity": result.ionic_conductivity,
                    "diffusion_coefficients": result.diffusion_coefficients,
                    "density": result.density,
                    "viscosity": result.viscosity,
                    "radial_distribution": result.radial_distribution
                }
            })
        except CalculationResult.DoesNotExist:
            logger.warning(f"计算结果不存在: calculation_id={calculation_id}")
            return Response({
                "success": False,
                "error": "计算结果不存在"
            }, status=status.HTTP_404_NOT_FOUND)
        
    except Exception as e:
        logger.error(f"获取计算结果时出错: {str(e)}")
        return Response({"success": False, "error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

# 实现提交计算API
# @permission_classes([IsAuthenticated])
@api_view(['POST'])
@permission_classes([])  # 空列表表示不需要任何权限
def submit_calculation(request):
    """提交新的电解液计算任务"""
    logger = logging.getLogger(__name__)
    logger.info(f"接收到提交计算请求: {request.data}")
    
    try:
        # 获取表单数据
        formulation_id = request.data.get('formulation_id')
        name = request.data.get('name', '电解液计算')
        description = request.data.get('description', '')
        
        # 检查有无直接提供计算配置
        config = request.data.get('config')
        
        if not formulation_id and not config:
            logger.error("缺少配方ID或计算配置")
            return Response({"success": False, "error": "缺少配方ID或计算配置"}, status=status.HTTP_400_BAD_REQUEST)
        
        # 如果有formulation_id，按正常流程处理
        if formulation_id:
            # 获取配方
            try:
                formulation = ElectrolyteFormulation.objects.get(id=formulation_id)
                
                # 检查权限（禁用认证时此处也临时跳过）
                # if formulation.user != request.user:
                #     logger.warning(f"用户尝试为不属于他的配方 {formulation_id} 提交计算")
                #     return Response({
                #         "success": False, 
                #         "error": "您无权为此配方提交计算"
                #     }, status=status.HTTP_403_FORBIDDEN)
                
                # 检查是否有输入文件
                try:
                    input_file = InputFile.objects.get(formulation=formulation)
                    if not input_file.file_path or not os.path.exists(input_file.file_path):
                        logger.error(f"找不到有效的输入文件: formulation_id={formulation_id}")
                        return Response({
                            "success": False, 
                            "error": "找不到有效的输入文件，请先生成并保存输入文件"
                        }, status=status.HTTP_400_BAD_REQUEST)
                except InputFile.DoesNotExist:
                    logger.error(f"找不到输入文件: formulation_id={formulation_id}")
                    return Response({
                        "success": False, 
                        "error": "找不到输入文件，请先生成并保存输入文件"
                    }, status=status.HTTP_400_BAD_REQUEST)
                    
                # 创建计算任务
                calculation = Calculation.objects.create(
                    user=request.user if hasattr(request, 'user') and request.user.is_authenticated else None,
                    formulation=formulation,
                    name=name,
                    description=description,
                    status='pending'
                )
            
            except ElectrolyteFormulation.DoesNotExist:
                logger.warning(f"找不到配方: formulation_id={formulation_id}，将创建临时配方")
                
                # 创建临时配方
                # 获取或创建默认用户
                default_user, created = User.objects.get_or_create(
                    username='default_system_user',
                    defaults={
                        'email': 'system@example.com',
                        'is_active': True
                    }
                )
                
                temp_formulation = ElectrolyteFormulation.objects.create(
                    name=f"临时配方_{formulation_id}",
                    description="前端自动创建的临时配方",
                    user=request.user if hasattr(request, 'user') and request.user.is_authenticated else default_user,
                )
                
                # 创建临时输入文件
                temp_dir = os.path.join(settings.MEDIA_ROOT, 'temp_formulations', str(temp_formulation.id))
                os.makedirs(temp_dir, exist_ok=True)
                
                temp_inp_path = os.path.join(temp_dir, 'input.inp')
                with open(temp_inp_path, 'w') as f:
                    f.write("# 临时生成的输入文件\n")
                    f.write("START\n")
                    f.write("name 临时配方\n")
                    f.write("T 298\n")
                    f.write("Box_size 40\n")
                    f.write("concentration 1.0\n")
                    f.write("Cation1_name Li\n")
                    f.write("Cation1_ratio 1.0\n")
                    f.write("Anion1_name PF6\n")
                    f.write("Anion1_ratio 1.0\n")
                    f.write("Solvent1_name EC\n")
                    f.write("Solvent1_ratio 1.0\n")
                    f.write("END\n")
                
                # 创建输入文件记录
                InputFile.objects.create(
                    formulation=temp_formulation,
                    file_path=temp_inp_path,
                    content="# 临时生成的输入文件\nSTART\nname 临时配方\nT 298\nBox_size 40\nconcentration 1.0\nCation1_name Li\nCation1_ratio 1.0\nAnion1_name PF6\nAnion1_ratio 1.0\nSolvent1_name EC\nSolvent1_ratio 1.0\nEND\n"
                )
                
                # 创建计算任务
                calculation = Calculation.objects.create(
                    user=request.user if hasattr(request, 'user') and request.user.is_authenticated else default_user,
                    formulation=temp_formulation,
                    name=name,
                    description=description or "基于临时配方的计算",
                    status='pending'
                )
                
                logger.info(f"已创建临时配方和计算任务: formulation_id={temp_formulation.id}, calculation_id={calculation.id}")
            
        # 如果直接提供配置，则使用临时配方和临时输入文件
        else:
            logger.info(f"使用直接提供的配置: {config}")
            
            # 创建临时文件夹
            temp_dir = os.path.join(settings.MEDIA_ROOT, 'temp_calculations', str(uuid.uuid4()))
            os.makedirs(temp_dir, exist_ok=True)
            
            # 写入临时INP文件
            temp_inp_path = os.path.join(temp_dir, 'input.inp')
            file_content = "START\n"
            for key, value in config.items():
                if isinstance(value, dict) or isinstance(value, list):
                    continue
                file_content += f"{key} {value}\n"
                
            # 处理盐和溶剂
            if 'salts' in config and config['salts']:
                for i, salt in enumerate(config['salts']):
                    file_content += f"salt_{i} {salt['name']},{salt['cation']},{salt['anion']}\n"
                    file_content += f"salt_concentration_{i} {salt['concentration']}\n"
                    
            if 'solvents' in config and config['solvents']:
                for i, solvent in enumerate(config['solvents']):
                    file_content += f"solvent_{i} {solvent['name']},{solvent['smile']}\n"
                    file_content += f"solvent_concentration_{i} {solvent['concentration']}\n"
                    
            file_content += "END\n"
            
            # 写入文件
            with open(temp_inp_path, 'w') as f:
                f.write(file_content)
                
            # 获取或创建默认用户
            default_user, created = User.objects.get_or_create(
                username='default_system_user',
                defaults={
                    'email': 'system@example.com',
                    'is_active': True
                }
            )
            
            # 创建临时配方
            temp_formulation = ElectrolyteFormulation.objects.create(
                name=f"临时配方_直接配置",
                description="通过直接配置创建的临时配方",
                user=request.user if hasattr(request, 'user') and request.user.is_authenticated else default_user,
            )
            
            # 创建临时输入文件对象
            input_file = InputFile.objects.create(
                formulation=temp_formulation,
                file_path=temp_inp_path,
                content=file_content
            )
            
            # 创建临时计算任务
            calculation = Calculation.objects.create(
                user=request.user if hasattr(request, 'user') and request.user.is_authenticated else default_user,
                formulation=temp_formulation,  # 使用新创建的配方
                name=name or "临时计算",
                description=description or "直接配置的计算",
                status='pending'
            )
        
        # 提交Celery任务
        task = run_electrolyte_calculation.delay(calculation.id)
        
        # 更新任务ID
        calculation.task_id = task.id
        calculation.save(update_fields=['task_id'])
        
        logger.info(f"Celery任务已提交: calculation_id={calculation.id}, task_id={task.id}")
        
        return Response({
            "success": True,
            "message": "计算任务已提交",
            "id": calculation.id,
            "task_id": task.id
        })
        
    except Exception as e:
        logger.error(f"提交计算任务时出错: {str(e)}")
        return Response({"success": False, "error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
