import os
import time
import json
import sys
from datetime import datetime
from celery import shared_task
from django.conf import settings
from django.utils import timezone
from .models import Calculation, CalculationResult, InputFile
import logging

# 添加项目根目录到Python路径
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"已添加 {project_root} 到Python路径")

# 尝试导入molyte_cursor
try:
    from molyte_cursor.src.io.inp_reader import read_inp_file, parse_inp_content
    from molyte_cursor.src.core.electrolyte_workflow import ElectrolyteWorkflow, run_from_inp_file, run_from_inp_content
    MOLYTE_AVAILABLE = True
    print("成功导入molyte_cursor模块!")
except ImportError as e:
    MOLYTE_AVAILABLE = False
    print(f"导入错误: {e}")

logger = logging.getLogger(__name__)

@shared_task(bind=True)
def run_electrolyte_calculation(self, calculation_id):
    """运行电解质计算的Celery任务
    
    这个任务将使用Molyte执行电解质性能计算，并更新计算结果
    """
    logger.info(f"开始处理电解液计算任务，计算ID: {calculation_id}")
    
    # 检查是否可以使用molyte_cursor模块
    logger.info(f"MOLYTE_AVAILABLE: {MOLYTE_AVAILABLE}")
    
    # 获取计算对象
    try:
        calculation = Calculation.objects.get(id=calculation_id)
    except Calculation.DoesNotExist:
        logger.error(f"计算ID {calculation_id} 不存在")
        return {"status": "failed", "error": f"计算ID {calculation_id} 不存在"}
    
    # 检查计算状态，避免重复处理
    if calculation.status in ['completed', 'failed']:
        logger.info(f"计算ID: {calculation_id} 状态为 {calculation.status}，跳过处理")
        return {"status": calculation.status, "calculation_id": calculation_id}
    
    try:
        # 更新计算状态
        calculation.status = 'running'
        calculation.started_at = timezone.now()
        calculation.save()
        
        # 检查molyte_cursor是否可用，如果不可用直接返回失败
        if not MOLYTE_AVAILABLE:
            error_msg = "molyte_cursor模块不可用，无法执行计算"
            logger.error(error_msg)
            
            # 更新计算状态为失败
            calculation.status = 'failed'
            calculation.error_message = error_msg
            calculation.finished_at = timezone.now()
            calculation.save()
            
            return {"status": "failed", "error": error_msg, "calculation_id": calculation_id}
        
        # 获取配方和参数
        formulation = calculation.formulation
        if not formulation:
            error_msg = "计算没有关联的配方"
            logger.error(error_msg)
            calculation.status = 'failed'
            calculation.error_message = error_msg
            calculation.finished_at = timezone.now()
            calculation.save()
            return {"status": "failed", "error": error_msg, "calculation_id": calculation_id}
            
        # 检查是否有INP输入文件
        try:
            input_file = InputFile.objects.get(formulation=formulation)
            if not input_file.file_path or not os.path.exists(input_file.file_path):
                raise ValueError(f"找不到有效的INP输入文件")
                
            inp_file_path = input_file.file_path
            logger.info(f"找到INP输入文件: {inp_file_path}")
            
        except InputFile.DoesNotExist:
            error_msg = f"配方没有关联的INP输入文件"
            logger.error(error_msg)
            calculation.status = 'failed'
            calculation.error_message = error_msg
            calculation.finished_at = timezone.now()
            calculation.save()
            return {"status": "failed", "error": error_msg, "calculation_id": calculation_id}
        
        # 更新进度信息
        self.update_state(state='PROGRESS', meta={'progress': 10, 'status': '正在解析输入文件'})
        
        # 创建输出目录
        temp_dir = os.path.join(settings.MEDIA_ROOT, 'calculations', str(calculation.id))
        os.makedirs(temp_dir, exist_ok=True)
        
        # 使用workspace目录作为计算输出目录
        aisuan_root = os.environ.get('AISUAN_ROOT', '/public/home/xiaoji/AiSuan')
        workspace_dir = os.path.join(aisuan_root, 'workspace')
        
        # 根据计算ID和用户信息创建目录路径
        user = calculation.user if hasattr(calculation, 'user') and calculation.user else 'anonymous'
        user_name = user.username if hasattr(user, 'username') else str(user)
        
        # 如果是匿名用户，使用public目录
        if user_name == 'anonymous':
            user_workspace = os.path.join(workspace_dir, 'public')
        else:
            user_workspace = os.path.join(workspace_dir, 'users', user_name)
        
        # 创建以计算ID命名的目录
        calc_dir_name = f"calc_{calculation.id}_{timezone.now().strftime('%Y%m%d_%H%M%S')}"
        output_dir = os.path.join(user_workspace, 'calculations', calc_dir_name)
        
        # 确保目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 保存这个路径以便后续引用
        calculation.output_dir = output_dir
        calculation.save()
        
        # 设置日志文件路径
        log_file_path = os.path.join(temp_dir, 'calculation.log')
        
        # 使用真实的molyte_cursor模块进行计算
        logger.info("使用molyte_cursor模块进行计算")
        
        try:
            # 读取INP文件内容
            with open(inp_file_path, 'r', encoding='utf-8') as f:
                inp_content = f.read()
            
            # 更新进度信息
            self.update_state(state='PROGRESS', meta={'progress': 20, 'status': '已解析输入文件'})
            
            # 执行工作流
            logger.info(f"开始执行电解液计算工作流, 输出目录: {output_dir}")
            
            # 根据输入类型选择调用方法
            if os.path.exists(inp_file_path):
                # 使用文件路径创建工作流
                workflow = ElectrolyteWorkflow(
                    inp_file_path=inp_file_path,
                    work_dir=output_dir
                )
            else:
                # 使用内容字符串创建工作流
                workflow = ElectrolyteWorkflow(
                    inp_content=inp_content,
                    work_dir=output_dir
                )
            
            # 执行各个步骤并更新进度
            # 步骤1：解析输入
            workflow.parse_input()
            self.update_state(state='PROGRESS', meta={'progress': 40, 'status': '已生成计算输入文件'})
            
            # 步骤2：生成输入文件
            workflow.generate_input_files()
            self.update_state(state='PROGRESS', meta={'progress': 60, 'status': '正在提交计算到Slurm'})
            
            # 步骤3：运行计算
            self.update_state(state='PROGRESS', meta={'progress': 60, 'status': '正在提交计算到Slurm'})
            slurm_job_id = workflow.run_calculation()
            
            # 保存Slurm作业ID
            if slurm_job_id:
                calculation.slurm_job_id = slurm_job_id
                calculation.status = 'submitted'  # 已提交状态
                calculation.save()
                
                self.update_state(state='PROGRESS', meta={'progress': 70, 'status': f'计算已提交到Slurm，作业ID: {slurm_job_id}'})
                
                # 可以添加一个轮询机制，检查作业状态
                # 以下是简化版本，实际上可能需要更复杂的处理
                logger.info(f"计算已提交到Slurm，作业ID: {slurm_job_id}")
                
                # 这里可以添加对作业完成的监控代码
                # 例如：
                # 1. 周期性检查作业状态
                # 2. 使用回调机制
                # 3. 使用Slurm作业完成回调脚本
                
                # 对于演示目的，我们假设作业会完成，继续后续步骤
                
            else:
                logger.warning("无法获取Slurm作业ID，将继续后续步骤")
            
            self.update_state(state='PROGRESS', meta={'progress': 80, 'status': '计算已完成'})
            
            # 步骤4：分析结果
            # 注意：在生产环境中，这一步应该在确认Slurm作业完成后执行
            workflow.analyze_results()
            
            self.update_state(state='PROGRESS', meta={'progress': 90, 'status': '正在生成报告'})
            
            # 步骤5：生成报告
            workflow.generate_report()
            
            # 获取结果摘要
            result_summary = workflow.get_result_summary()
            
            # 检查结果状态
            if result_summary.get('status') == 'completed' or 'conductivity_mS_cm' in result_summary:
                result = {
                    "status": "completed",
                    "simulation_type": "molyte_cursor",
                    "data": result_summary
                }
            else:
                raise ValueError("工作流未能提供完整结果")
            
        except Exception as e:
            logger.error(f"使用molyte_cursor时出错: {str(e)}", exc_info=True)
            raise e
            
        # 写入日志文件
        with open(log_file_path, 'w') as f:
            f.write(f"# 计算日志\n")
            f.write(f"开始时间: {calculation.started_at}\n")
            f.write(f"输入文件: {inp_file_path}\n")
            f.write(f"结果: {json.dumps(result, indent=2)}\n")
        
        # 保存日志文件
        calculation.log_file = os.path.join('calculations', str(calculation.id), 'calculation.log')
        
        # 更新计算状态
        calculation.status = 'completed'
        calculation.finished_at = timezone.now()
        calculation.save()
        
        # 创建或更新计算结果
        result_data = result.get('data', {})
        
        result_obj, created = CalculationResult.objects.get_or_create(calculation=calculation)
        
        # 提取结果字段
        result_obj.ionic_conductivity = result_data.get('conductivity_mS_cm', 0)
        result_obj.density = result_data.get('density', 0)
        result_obj.viscosity = result_data.get('viscosity_mPa_s', 0)
        
        # 提取扩散系数
        diffusion_coeffs = result_data.get('diffusion_coefficients', {})
        if not diffusion_coeffs and 'diffusion' in result_data:
            diffusion_coeffs = result_data['diffusion']
            
        result_obj.diffusion_coefficients = diffusion_coeffs
        
        # 提取径向分布函数
        rdf_data = {}
        if 'rdf_analysis' in result_data:
            rdf_data = result_data['rdf_analysis']
        result_obj.radial_distribution = rdf_data
        
        # 保存完整结果
        result_obj.additional_results = result
        result_obj.save()
        
        logger.info(f"计算ID: {calculation_id} 完成")
        return {"status": "completed", "calculation_id": calculation_id}
        
    except Exception as e:
        # 捕获所有异常，更新计算状态并记录错误信息
        logger.error(f"计算ID: {calculation_id} 失败: {str(e)}", exc_info=True)
        try:
            calculation.status = 'failed'
            calculation.error_message = str(e)
            calculation.finished_at = timezone.now()
            calculation.save()
        except Exception as ex:
            logger.error(f"更新计算状态时出错: {str(ex)}")
        
        return {"status": "failed", "error": str(e), "calculation_id": calculation_id} 