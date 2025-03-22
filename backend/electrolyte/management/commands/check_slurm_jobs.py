"""
Slurm作业状态检查命令

检查已提交的Slurm作业状态并更新数据库
"""
import os
import sys
import logging
from django.core.management.base import BaseCommand
from electrolyte.models import Calculation
from django.utils import timezone

# 添加项目根目录到Python路径
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# 尝试导入Slurm监控工具
try:
    from molyte_cursor.src.utils.slurm_monitor import check_job_status, get_job_info, get_job_log
    MONITOR_AVAILABLE = True
except ImportError as e:
    MONITOR_AVAILABLE = False
    print(f"导入Slurm监控工具失败: {e}")

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = '检查已提交的Slurm作业状态并更新数据库'
    
    def add_arguments(self, parser):
        parser.add_argument(
            '--job-id',
            help='指定要检查的Slurm作业ID',
        )
        parser.add_argument(
            '--calculation-id',
            help='指定要检查的计算ID',
        )
        parser.add_argument(
            '--all',
            action='store_true',
            help='检查所有已提交但未完成的作业',
        )
    
    def handle(self, *args, **options):
        if not MONITOR_AVAILABLE:
            self.stderr.write(self.style.ERROR('Slurm监控工具不可用，无法执行命令'))
            return
        
        # 检查参数
        job_id = options.get('job_id')
        calculation_id = options.get('calculation_id')
        check_all = options.get('all')
        
        if job_id:
            # 根据Slurm作业ID查询
            calculations = Calculation.objects.filter(slurm_job_id=job_id)
            if not calculations.exists():
                self.stderr.write(self.style.WARNING(f'找不到Slurm作业ID为 {job_id} 的计算'))
                return
            
            for calculation in calculations:
                self._check_and_update_job(calculation)
                
        elif calculation_id:
            # 根据计算ID查询
            try:
                calculation = Calculation.objects.get(id=calculation_id)
                if not calculation.slurm_job_id:
                    self.stderr.write(self.style.WARNING(f'计算ID为 {calculation_id} 的计算没有关联的Slurm作业ID'))
                    return
                
                self._check_and_update_job(calculation)
                
            except Calculation.DoesNotExist:
                self.stderr.write(self.style.ERROR(f'找不到计算ID为 {calculation_id} 的计算'))
                return
                
        elif check_all:
            # 检查所有已提交但未完成的作业
            pending_calculations = Calculation.objects.filter(
                status__in=['submitted', 'queued', 'running'],
                slurm_job_id__isnull=False
            )
            
            if not pending_calculations.exists():
                self.stdout.write(self.style.SUCCESS('没有待处理的Slurm作业'))
                return
            
            self.stdout.write(f'找到 {pending_calculations.count()} 个待处理的Slurm作业')
            
            for calculation in pending_calculations:
                self._check_and_update_job(calculation)
                
        else:
            self.stderr.write(self.style.ERROR('请指定要检查的作业（--job-id, --calculation-id 或 --all）'))
    
    def _check_and_update_job(self, calculation):
        """检查作业状态并更新计算记录"""
        job_id = calculation.slurm_job_id
        
        self.stdout.write(f'检查计算ID {calculation.id} 的Slurm作业 {job_id}')
        
        # 检查作业状态
        status = check_job_status(job_id)
        
        self.stdout.write(f'Slurm作业 {job_id} 状态: {status}')
        
        # 记录上次检查状态
        previous_status = calculation.status
        
        # 根据Slurm状态更新计算状态
        if status == 'PENDING':
            # 作业仍在排队
            calculation.status = 'queued'
            self.stdout.write(self.style.WARNING(f'作业 {job_id} 仍在排队'))
            
        elif status == 'RUNNING':
            # 作业正在运行
            calculation.status = 'running'
            self.stdout.write(self.style.WARNING(f'作业 {job_id} 正在运行'))
            
        elif status == 'COMPLETED':
            # 作业已成功完成
            calculation.status = 'completed'
            calculation.finished_at = timezone.now()
            self.stdout.write(self.style.SUCCESS(f'作业 {job_id} 已成功完成'))
            
            try:
                # 获取作业日志
                output_log = get_job_log(job_id, 'out')
                error_log = get_job_log(job_id, 'err')
                
                # 获取作业完成日志
                job_completed_log = os.path.join(calculation.output_dir, 'job_completed.log')
                completed_log_content = ""
                if os.path.exists(job_completed_log):
                    with open(job_completed_log, 'r', encoding='utf-8', errors='replace') as f:
                        completed_log_content = f.read()
                
                # 保存所有日志
                combined_log_path = os.path.join(calculation.output_dir, 'slurm_combined.log')
                with open(combined_log_path, 'w', encoding='utf-8') as f:
                    f.write(f"# Slurm作业 {job_id} 综合日志\n\n")
                    f.write("## 作业完成信息\n")
                    f.write(completed_log_content)
                    f.write("\n\n## 标准输出\n")
                    f.write(output_log)
                    f.write("\n\n## 标准错误\n")
                    f.write(error_log)
                
                self.stdout.write(f'已保存综合日志到 {combined_log_path}')
                
                # 检查是否存在结果文件目录
                results_dir = os.path.join(calculation.output_dir, 'results')
                if os.path.exists(results_dir) and os.path.isdir(results_dir):
                    self.stdout.write(f'结果文件目录存在: {results_dir}')
                    
                    # 检查log.lammps文件
                    lammps_log = os.path.join(results_dir, 'log.lammps')
                    if os.path.exists(lammps_log):
                        self.stdout.write(self.style.SUCCESS(f'找到LAMMPS日志文件: {lammps_log}'))
                    else:
                        self.stdout.write(self.style.WARNING(f'未找到LAMMPS日志文件'))
                        
                        # 尝试从工作目录复制
                        work_dir_log = None
                        for root, dirs, files in os.walk(calculation.output_dir):
                            if 'log.lammps' in files:
                                work_dir_log = os.path.join(root, 'log.lammps')
                                break
                        
                        if work_dir_log:
                            import shutil
                            shutil.copy(work_dir_log, lammps_log)
                            self.stdout.write(self.style.SUCCESS(f'已复制LAMMPS日志文件到结果目录'))
                            
                # 导入必要模块以进行结果分析
                try:
                    from molyte_cursor.src.core.electrolyte_workflow import ElectrolyteWorkflow
                    
                    # 创建工作流对象并加载已有结果
                    workflow = ElectrolyteWorkflow(output_dir=calculation.output_dir)
                    
                    # 分析结果
                    try:
                        self.stdout.write('开始分析计算结果...')
                        workflow.analyze_results()
                        workflow.generate_report()
                        result_summary = workflow.get_result_summary()
                        
                        # 更新计算结果
                        from electrolyte.models import CalculationResult
                        result_obj, created = CalculationResult.objects.get_or_create(calculation=calculation)
                        
                        # 提取结果字段
                        result_data = result_summary
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
                        result_obj.additional_results = {"status": "completed", "data": result_summary}
                        result_obj.save()
                        
                        self.stdout.write(self.style.SUCCESS(f'已成功分析并保存计算结果'))
                        
                    except Exception as e:
                        self.stdout.write(self.style.ERROR(f'分析结果时出错: {str(e)}'))
                        calculation.error_message = f"分析结果时出错: {str(e)}"
                except ImportError as e:
                    self.stdout.write(self.style.ERROR(f'导入工作流模块失败: {str(e)}'))
                
            except Exception as e:
                self.stdout.write(self.style.ERROR(f'处理作业日志时出错: {str(e)}'))
                calculation.error_message = f"处理作业日志时出错: {str(e)}"
            
        elif status == 'FAILED':
            # 作业失败
            calculation.status = 'failed'
            calculation.finished_at = timezone.now()
            calculation.error_message = f"Slurm作业 {job_id} 失败"
            self.stdout.write(self.style.ERROR(f'作业 {job_id} 失败'))
            
            # 获取错误日志
            error_log = get_job_log(job_id, 'err')
            if error_log:
                calculation.error_message += f"\n\n{error_log}"
                
            # 检查作业失败原因
            job_info = get_job_info(job_id)
            if 'JobID' in job_info and 'State' in job_info:
                failure_state = job_info['State']
                calculation.error_message += f"\n\nSlurm状态: {failure_state}"
                self.stdout.write(self.style.ERROR(f'详细失败状态: {failure_state}'))
                
            # 检查是否有作业完成日志
            job_completed_log = os.path.join(calculation.output_dir, 'job_completed.log')
            if os.path.exists(job_completed_log):
                with open(job_completed_log, 'r', encoding='utf-8', errors='replace') as f:
                    completed_log_content = f.read()
                    calculation.error_message += f"\n\n作业完成日志:\n{completed_log_content}"
            
        elif status == 'CANCELLED':
            # 作业被取消
            calculation.status = 'cancelled'
            calculation.finished_at = timezone.now()
            self.stdout.write(self.style.WARNING(f'作业 {job_id} 被取消'))
            
            # 获取取消原因
            job_info = get_job_info(job_id)
            if 'JobID' in job_info and 'State' in job_info:
                cancel_state = job_info['State']
                calculation.error_message = f"作业被取消，状态: {cancel_state}"
                
        else:
            # 未知状态
            self.stdout.write(self.style.WARNING(f'作业 {job_id} 状态未知: {status}'))
            
            # 检查作业是否存在
            if status == "UNKNOWN":
                # 检查作业是否已运行超过7天
                if calculation.started_at and (timezone.now() - calculation.started_at).days > 7:
                    self.stdout.write(self.style.WARNING(f'作业 {job_id} 已提交超过7天且状态未知，标记为失败'))
                    calculation.status = 'failed'
                    calculation.finished_at = timezone.now()
                    calculation.error_message = f"作业提交后超过7天无法获取状态，可能已丢失"
        
        # 检查状态是否发生变化
        if previous_status != calculation.status:
            self.stdout.write(self.style.SUCCESS(f'计算状态从 {previous_status} 更新为 {calculation.status}'))
            
        # 保存更新后的计算记录
        calculation.save()
        self.stdout.write(f'已更新计算ID {calculation.id} 的状态为 {calculation.status}') 