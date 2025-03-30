import os
import time
import json
import sys
import traceback  # 添加traceback模块导入
from datetime import datetime
from celery import shared_task
from django.conf import settings
from django.utils import timezone
from .models import Calculation, CalculationResult, InputFile, ElectrolyteFormulation, FormulationComponent, SimulationParameters
import logging

# 添加项目根目录到Python路径
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"已添加 {project_root} 到Python路径")

# 创建日志记录器
logger = logging.getLogger(__name__)

# 确保可以找到molyte_cursor模块
molyte_paths = [
    os.path.join(project_root, 'molyte_cursor'),
    '/public/home/xiaoji/molyte_v1/molyte_cursor',
    '/public/home/xiaoji/AiSuan/molyte_cursor',
]

for path in molyte_paths:
    if os.path.exists(path) and path not in sys.path:
        sys.path.insert(0, path)
        logger.info(f"添加molyte_cursor路径: {path}")
        print(f"添加molyte_cursor路径: {path}")

# 尝试导入molyte_cursor
try:
    from src.io.inp_reader import read_inp_file, parse_inp_content
    # 尝试导入新的工作流
    try:
        from src.core.molyte_workflow import MolyteWorkflow, run_from_inp_file, run_from_inp_content
        # 注释掉动态添加方法的代码，因为MolyteWorkflow现在已经有这个方法了
        # if not hasattr(MolyteWorkflow, 'from_formulation'):
        #     def from_formulation(cls, formulation_data=None, work_dir=None, output_dir=None,
        #                         user_id=None, project_id=None, formulation_id=None,
        #                         task_id=None, config=None, **kwargs):
        #         """从配方数据创建工作流实例"""
        #         # 创建实例
        #         instance = cls(
        #             inp_file_path=None,
        #             inp_content=None,
        #             work_dir=work_dir,
        #             output_dir=output_dir,
        #             user_id=user_id,
        #             project_id=project_id,
        #             formulation_id=formulation_id,
        #             task_id=task_id,
        #             formulation_data=formulation_data,
        #             config=config
        #         )
        #         return instance
        #     
        #     # 添加为类方法
        #     MolyteWorkflow.from_formulation = classmethod(from_formulation)
        #     print("添加了from_formulation类方法")
        MOLYTE_WORKFLOW_AVAILABLE = True
        print("成功导入MolyteWorkflow模块!")
    except ImportError as e:
        MOLYTE_WORKFLOW_AVAILABLE = False
        print(f"导入MolyteWorkflow错误: {e}")
        # 回退到原始工作流
    MOLYTE_AVAILABLE = True
    print("成功导入molyte_cursor模块!")
except ImportError as e:
    MOLYTE_WORKFLOW_AVAILABLE = False
    MOLYTE_AVAILABLE = False
    print(f"导入错误: {e}")
    # 尝试直接从绝对路径导入
    try:
        sys.path.insert(0, '/public/home/xiaoji/AiSuan')
        from molyte_cursor.src.core.molyte_workflow import MolyteWorkflow, run_from_inp_file
        from molyte_cursor.src.io.inp_reader import read_inp_file, parse_inp_content
        MOLYTE_WORKFLOW_AVAILABLE = True
        MOLYTE_AVAILABLE = True
        print("使用绝对路径成功导入molyte_cursor模块!")
    except ImportError as e:
        print(f"使用绝对路径导入失败: {e}")

if not MOLYTE_WORKFLOW_AVAILABLE:
    class MockMolyteWorkflow:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
        
        def run(self):
            return {"status": "failed", "error": "MolyteWorkflow模块不可用"}
    
    MolyteWorkflow = MockMolyteWorkflow
    print("使用Mock MolyteWorkflow")

@shared_task(bind=True)
def run_electrolyte_calculation(self, calculation_id, auto_submit=True):
    """运行电解液计算任务
    
    Args:
        calculation_id: 计算ID
        auto_submit: 是否自动提交SLURM作业
    """
    from django.apps import apps
    import os
    import time
    
    Calculation = apps.get_model('electrolyte', 'Calculation')
    logger.info(f"开始执行Celery任务: calculation_id={calculation_id}, auto_submit={auto_submit}")
    
    try:
        # 获取计算记录
        calculation = Calculation.objects.get(id=calculation_id)
        
        # 初始化工作目录和输出目录
        base_dir = os.environ.get('MOLYTE_OUTPUT_DIR', '/public/home/xiaoji/AiSuan/molyte_cursor')
        timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        work_dir = os.path.join(base_dir, "workspace", "users", "default", f"task_{timestamp}")
        output_dir = os.path.join(work_dir, "output")
        
        logger.info(f"计算任务详情: work_dir={work_dir}, output_dir={output_dir}")
        
        # 创建工作流实例
        try:
            logger.info(f"创建工作流实例: calculation_id={calculation_id}")
            workflow = _create_workflow_from_id(
                calculation_id=calculation_id,
                work_dir=work_dir,
                output_dir=output_dir
            )
            
            if not workflow:
                raise ValueError("工作流实例创建失败")
            
            logger.info(f"工作流实例创建成功: {workflow}")
        except Exception as e:
            error_msg = f"创建工作流实例失败: {str(e)}"
            logger.error(error_msg, exc_info=True)
            
            # 更新计算记录状态为失败
            calculation.status = 'failed'
            calculation.error_message = error_msg[:1000]  # 限制错误消息长度
            calculation.save(update_fields=['status', 'error_message'])
            
            raise ValueError(error_msg)
            
        # 创建输入文件
        try:
            logger.info("开始生成LAMMPS输入文件")
            workflow.create_lammps_files()
            logger.info("LAMMPS输入文件生成成功")
        except Exception as e:
            error_msg = f"生成LAMMPS输入文件失败: {str(e)}"
            logger.error(error_msg, exc_info=True)
            
            # 更新计算记录状态为失败
            calculation.status = 'failed'
            calculation.error_message = error_msg[:1000]
            calculation.save(update_fields=['status', 'error_message'])
            
            raise ValueError(error_msg)
            
        # 如果设置了自动提交，则提交SLURM作业
        if auto_submit:
            try:
                logger.info("开始提交SLURM作业")
                # 定义回调函数，用于更新数据库中的作业ID
                def update_job_id(job_id):
                    logger.info(f"收到作业ID回调: job_id={job_id}")
                    calculation.slurm_job_id = job_id
                    calculation.save(update_fields=['slurm_job_id'])
                    logger.info(f"成功更新数据库中的SLURM作业ID: {job_id}")
                
                # 在提交作业时传递回调函数
                job_id = workflow.submit_job(callback=update_job_id)
                
                if job_id:
                    logger.info(f"SLURM作业提交成功: job_id={job_id}")
                    
                    # 确保数据库中的作业ID已更新（防止回调未执行）
                    if not calculation.slurm_job_id:
                        calculation.slurm_job_id = job_id
                        calculation.save(update_fields=['slurm_job_id'])
                        logger.info(f"手动更新数据库中的SLURM作业ID: {job_id}")
                    
                    return {
                        "calculation_id": calculation_id,
                        "status": "running",
                        "slurm_job_id": job_id,
                        "work_dir": work_dir,
                        "output_dir": output_dir
                    }
                else:
                    error_msg = "SLURM作业提交失败: 未获取到有效的作业ID"
                    logger.error(error_msg)
                    
                    # 更新计算记录状态为失败
                    calculation.status = 'failed'
                    calculation.error_message = error_msg
                    calculation.save(update_fields=['status', 'error_message'])
                    
                    raise ValueError(error_msg)
                    
            except Exception as e:
                error_msg = f"提交SLURM作业失败: {str(e)}"
                logger.error(error_msg, exc_info=True)
                
                # 更新计算记录状态为失败
                calculation.status = 'failed'
                calculation.error_message = error_msg[:1000]
                calculation.save(update_fields=['status', 'error_message'])
                
                raise ValueError(error_msg)
        else:
            logger.info("自动提交SLURM作业被禁用，需要手动提交")
            return {
                "calculation_id": calculation_id,
                "status": "prepared",
                "slurm_job_id": None,
                "work_dir": work_dir,
                "output_dir": output_dir
            }
    except Calculation.DoesNotExist:
        error_msg = f"计算ID {calculation_id} 不存在"
        logger.error(error_msg)
        raise ValueError(error_msg)
    except Exception as e:
        error_msg = f"执行计算任务时发生未预期的错误: {str(e)}"
        logger.error(error_msg, exc_info=True)
        
        try:
            calculation = Calculation.objects.get(id=calculation_id)
            calculation.status = 'failed'
            calculation.error_message = error_msg[:1000]
            calculation.save(update_fields=['status', 'error_message'])
        except:
            logger.error(f"无法更新计算ID {calculation_id} 的状态为失败", exc_info=True)
        
        raise ValueError(error_msg)

def _create_workflow_from_id(calculation_id, work_dir=None, output_dir=None):
    """从计算ID创建工作流实例"""
    # 导入必要的Django模型
    from django.apps import apps
    import os  # 在函数顶部导入os模块
    Calculation = apps.get_model('electrolyte', 'Calculation')
    InputFile = apps.get_model('electrolyte', 'InputFile')
    
    try:
        # 获取计算记录
        calculation = Calculation.objects.get(id=calculation_id)
        formulation = calculation.formulation
        
        if not formulation:
            error_msg = f"计算ID {calculation_id} 没有关联的配方"
            logger.error(error_msg)
            raise ValueError(error_msg)
        
        # 选择工作流创建方式 - 优先使用直接从配方创建的方式
        try:
            # 尝试直接从配方数据创建工作流
            logger.info(f"尝试直接从配方数据创建工作流: formulation_id={formulation.id}")
            
            # 获取SimulationParameters实例（如果存在）
            try:
                sim_params = SimulationParameters.objects.get(formulation_id=formulation.id)
                box_size = sim_params.box_size
                logger.info(f"【参数传递】从数据库读取的box_size: {box_size} ({type(box_size).__name__})")
            except SimulationParameters.DoesNotExist:
                logger.warning(f"未找到formulation_id={formulation.id}的SimulationParameters，使用默认值")
                box_size = 50.0  # 默认值
            
            # 确保box_size为float类型
            try:
                box_size = float(box_size)
                logger.info(f"【参数处理】转换后的box_size: {box_size} ({type(box_size).__name__})")
            except (ValueError, TypeError) as e:
                logger.error(f"【参数错误】box_size转换失败: {e}，使用默认值50.0")
                box_size = 50.0
            
            # 构建配方数据字典
            formulation_data = {
                "name": formulation.name,
                "description": formulation.description,
                "parameters": {
                    "temperature": sim_params.temperature if sim_params else 298.15,
                    "pressure": sim_params.pressure if sim_params else 1.0,
                    "box_size": box_size,
                    "cutoff": sim_params.cutoff if sim_params else 10.0,
                    "time_step": sim_params.time_step if sim_params else 1.0,
                    "equilibration_steps": sim_params.equilibration_steps if sim_params else 10000,
                    "production_steps": sim_params.production_steps if sim_params else 100000,
                },
                "cations": [],
                "anions": [],
                "solvents": [],
            }
            
            # 准备独立的config参数，包含auto_submit设置
            config = {
                "auto_submit": True,  # 启用自动提交
                "cpus": 64,           # 设置CPU数量
                "max_hours": 72       # 设置最大运行时间
            }
            
            logger.info(f"【参数检查】独立config设置: {config}")
            
            # 获取组件
            components = FormulationComponent.objects.filter(formulation=formulation)
            for component in components:
                if component.component_type == 'salt' and component.salt:
                    salt = component.salt
                    # 添加阳离子
                    formulation_data["cations"].append({
                        "name": salt.cation,
                        "number": component.number if component.number is not None else component.concentration,
                        "charge": 1,  # 默认电荷
                    })
                    # 添加阴离子
                    formulation_data["anions"].append({
                        "name": salt.anion,
                        "number": component.number if component.number is not None else component.concentration,
                        "charge": -1,  # 默认电荷
                    })
                elif component.component_type == 'solvent' and component.solvent:
                    solvent = component.solvent
                    # 添加溶剂
                    formulation_data["solvents"].append({
                        "name": solvent.name,
                        "smile": solvent.smile,
                        "number": component.number if component.number is not None else component.concentration,
                    })
            
            logger.info(f"已构建配方数据: {formulation_data}")
            
            # 检查构建的配方数据是否有有效的成分
            if not any([formulation_data['cations'], formulation_data['anions'], formulation_data['solvents']]):
                logger.warning("构建的配方数据没有有效的成分，可能会导致后续处理出错")
            
            try:
                # 使用配方数据创建工作流
                workflow = MolyteWorkflow.from_formulation(
                    formulation_data=formulation_data,
                    work_dir=work_dir,
                    output_dir=output_dir,
                    user_id=calculation.user.id if calculation.user else None,
                    project_id=None,
                    formulation_id=formulation.id,
                    task_id=run_electrolyte_calculation.request.id if hasattr(run_electrolyte_calculation, 'request') else None,
                    config=config  # 作为独立参数传递config
                )
                logger.info(f"成功直接从配方创建工作流，跳过读取INP文件步骤")
                
                # 更新计算任务的状态
                calculation.status = 'running'
                calculation.started_at = timezone.now()
                calculation.save(update_fields=['status', 'started_at'])
                
                return workflow
            except Exception as e:
                logger.error(f"创建MolyteWorkflow.from_formulation实例失败: {str(e)}")
                raise
            
        except Exception as e:
            # 直接创建失败，回退到使用INP文件
            logger.warning(f"直接从配方创建工作流失败: {str(e)}，尝试使用INP文件方式")
            
            # 获取输入文件
            try:
                input_file = InputFile.objects.get(formulation=formulation)
            except InputFile.DoesNotExist:
                # 如果没有输入文件，尝试生成一个
                logger.warning("未找到输入文件，尝试生成")
                
                # 导入必要的模块
                from .views import generate_lammps_inp_file
                from django.http import HttpRequest
                import os  # 确保在此处导入os模块
                
                # 创建一个基本的请求对象
                request = HttpRequest()
                request.user = calculation.user or None
                
                # 生成INP文件内容
                inp_content = generate_lammps_inp_file(formulation, request)
                
                if not inp_content:
                    error_msg = "无法生成INP文件内容"
                    logger.error(error_msg)
                    raise ValueError(error_msg)
                
                # 保存INP文件内容
                from django.conf import settings
                
                media_dir = settings.MEDIA_ROOT
                formulation_dir = os.path.join(media_dir, "formulations", f"formulation_{formulation.id}")
                os.makedirs(formulation_dir, exist_ok=True)
                
                file_path = os.path.join(formulation_dir, "电解液配方.inp")
                with open(file_path, "w") as f:
                    f.write(inp_content)
                
                # 创建InputFile记录
                from .models import InputFile
                relative_path = os.path.join("formulations", f"formulation_{formulation.id}", "电解液配方.inp")
                input_file = InputFile.objects.create(
                    formulation=formulation,
                    content=inp_content,
                    file_path=relative_path,
                    file_name="电解液配方.inp"
                )
                logger.info(f"已自动生成INP文件: {file_path}")
            
            # 构建完整的INP文件路径
            file_path = input_file.file_path
            # 检查是否已经是绝对路径
            if not os.path.isabs(file_path):
                # 添加MEDIA_ROOT前缀
                from django.conf import settings
                file_path = os.path.join(settings.MEDIA_ROOT, file_path)
                logger.info(f"添加MEDIA_ROOT前缀到文件路径: {file_path}")
                
            if not os.path.exists(file_path):
                error_msg = f"找不到有效的INP输入文件: {file_path} (数据库记录: {input_file.file_path})"
                logger.error(error_msg)
                raise ValueError(error_msg)
            
            # 使用文件路径而不是内容
            logger.info(f"成功获取输入文件: {file_path}")
            
            # 更新计算任务的状态
            calculation.status = 'running'
            calculation.started_at = timezone.now()
            calculation.save(update_fields=['status', 'started_at'])
            
            # 创建工作流，传递文件路径而不是内容
            # 准备config参数，包含更多设置
            config = {
                "auto_submit": True,  # 启用自动提交
                "cpus": 64,           # 设置CPU数量
                "max_hours": 72       # 设置最大运行时间
            }
            
            logger.info(f"【参数检查】使用INP文件方式时的config设置: {config}")
            
            # 创建一个默认的formulation_data，确保不会是None
            default_formulation_data = {
                "name": "电解液配方",
                "description": "电解液配方",
                "cations": [],
                "anions": [],
                "solvents": [],
                "parameters": {"box_size": 50.0}
            }
            
            return MolyteWorkflow(
                inp_file_path=file_path,  # 关键修改：使用绝对文件路径
                work_dir=work_dir,
                output_dir=output_dir,
                user_id=calculation.user.id if calculation.user else None,
                project_id=None,
                formulation_id=formulation.id,
                task_id=run_electrolyte_calculation.request.id if hasattr(run_electrolyte_calculation, 'request') else None,
                config=config,  # 完整的config参数
                formulation_data=default_formulation_data  # 添加默认的formulation_data
            )
        
    except Calculation.DoesNotExist:
        error_msg = f"计算ID {calculation_id} 不存在"
        logger.error(error_msg)
        raise ValueError(error_msg)
    except Exception as e:
        error_msg = f"创建工作流实例时发生未预期的错误: {str(e)}"
        logger.error(error_msg, exc_info=True)
        raise ValueError(error_msg)
