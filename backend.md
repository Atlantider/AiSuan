# 电解液计算模块后端设计方案

## 一、技术架构

### 后端框架
- **Web框架**: Django 4.2+ 和 Django REST Framework
- **任务队列**: Celery
- **消息代理/缓存**: Redis
- **API文档**: Swagger/OpenAPI
- **认证系统**: JWT
- **数据库**: PostgreSQL/MySQL

### 计算引擎集成
- **Molyte集成**: 通过Python API调用Molyte功能
- **LAMMPS集成**: 通过命令行调用执行分子动力学模拟
- **数据处理**: NumPy, pandas, scikit-learn
- **可视化生成**: Matplotlib, seaborn

## 二、系统架构图

```
前端 (React+Ant Design) <-> API层 (DRF) <-> 业务逻辑层 <-> Celery任务队列
                                                      |
                                                      v
                                              Molyte/LAMMPS计算引擎
                                                      |
                                                      v
                                                结果分析与存储
```

## 三、数据模型设计

### 核心数据模型

```python
# jobs/models.py
from django.db import models
from django.contrib.auth.models import User
import uuid

class Job(models.Model):
    """计算作业模型"""
    STATUS_CHOICES = (
        ('queued', '排队中'),
        ('running', '运行中'),
        ('completed', '已完成'),
        ('failed', '失败'),
    )
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='jobs')
    name = models.CharField(max_length=255)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='queued')
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    def __str__(self):
        return f"{self.name} ({self.status})"

class ElectrolyteFormula(models.Model):
    """电解液配方模型"""
    job = models.OneToOneField(Job, on_delete=models.CASCADE, related_name='formula')
    cations = models.JSONField(help_text="阳离子数据")
    anions = models.JSONField(help_text="阴离子数据")
    solvents = models.JSONField(help_text="溶剂数据")
    temperatures = models.JSONField(help_text="温度数据")
    concentrations = models.JSONField(help_text="浓度数据")
    calculation_types = models.JSONField(help_text="计算类型")
    simulation_parameters = models.JSONField(default=dict, help_text="模拟参数")
    
    def __str__(self):
        return f"Formula for {self.job.name}"

class SimulationResult(models.Model):
    """模拟结果模型"""
    job = models.OneToOneField(Job, on_delete=models.CASCADE, related_name='result')
    result_files = models.JSONField(help_text="结果文件路径")
    analysis_results = models.JSONField(help_text="分析结果数据")
    plots = models.JSONField(help_text="图表文件路径")
    summary = models.TextField(blank=True, help_text="结果摘要")
    
    def __str__(self):
        return f"Results for {self.job.name}"
```

### 辅助数据模型

```python
# electrolytes/models.py
from django.db import models

class CationOption(models.Model):
    """预定义阳离子选项"""
    label = models.CharField(max_length=100)
    value = models.CharField(max_length=50)
    charge = models.IntegerField()
    
    def __str__(self):
        return self.label

class AnionOption(models.Model):
    """预定义阴离子选项"""
    label = models.CharField(max_length=100)
    value = models.CharField(max_length=50)
    charge = models.IntegerField()
    smile = models.TextField(blank=True)
    
    def __str__(self):
        return self.label

class SolventOption(models.Model):
    """预定义溶剂选项"""
    label = models.CharField(max_length=100)
    value = models.CharField(max_length=50)
    smile = models.TextField()
    description = models.TextField(blank=True)
    
    def __str__(self):
        return self.label
```

## 四、API设计

### 电解液计算相关API

```
# 作业管理
POST /api/jobs/                       # 创建新作业
GET /api/jobs/                        # 获取作业列表
GET /api/jobs/{job_id}/               # 获取作业详情
DELETE /api/jobs/{job_id}/            # 删除作业

# 电解液计算
POST /api/electrolyte/calculate/      # 提交电解液计算
GET /api/electrolyte/jobs/{job_id}/status/   # 获取计算状态
GET /api/electrolyte/results/{job_id}/       # 获取计算结果

# 文件操作
POST /api/electrolyte/upload-template/       # 上传Excel模板数据
GET /api/electrolyte/download-template/      # 下载Excel模板
GET /api/electrolyte/download/{job_id}/{file_type}/  # 下载计算文件

# 选项数据
GET /api/electrolyte/options/cations/        # 获取阳离子选项
GET /api/electrolyte/options/anions/         # 获取阴离子选项
GET /api/electrolyte/options/solvents/       # 获取溶剂选项
```

### API实现示例

```python
# views.py
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response
from .models import Job, ElectrolyteFormula, SimulationResult
from .serializers import JobSerializer, ElectrolyteFormulaSerializer, SimulationResultSerializer
from .tasks import run_electrolyte_simulation

class JobViewSet(viewsets.ModelViewSet):
    """作业管理视图集"""
    queryset = Job.objects.all()
    serializer_class = JobSerializer
    
    def perform_create(self, serializer):
        serializer.save(user=self.request.user)

class ElectrolyteViewSet(viewsets.ViewSet):
    """电解液计算视图集"""
    
    @action(detail=False, methods=['post'])
    def calculate(self, request):
        """提交电解液计算"""
        job_serializer = JobSerializer(data={'name': request.data.get('name', '电解液计算')})
        if job_serializer.is_valid():
            job = job_serializer.save(user=request.user)
            
            formula_data = {
                'job': job.id,
                'cations': request.data.get('cations', []),
                'anions': request.data.get('anions', []),
                'solvents': request.data.get('solvents', []),
                'temperatures': request.data.get('temperatures', []),
                'concentrations': request.data.get('concentrations', []),
                'calculation_types': request.data.get('calculations', []),
                'simulation_parameters': request.data.get('parameters', {})
            }
            
            formula_serializer = ElectrolyteFormulaSerializer(data=formula_data)
            if formula_serializer.is_valid():
                formula_serializer.save()
                
                # 提交Celery任务
                run_electrolyte_simulation.delay(str(job.id))
                
                return Response({
                    'job_id': job.id,
                    'status': 'submitted',
                    'message': '计算任务已提交'
                }, status=status.HTTP_202_ACCEPTED)
            
            # 如果公式数据无效，删除作业并返回错误
            job.delete()
            return Response(formula_serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        
        return Response(job_serializer.errors, status=status.HTTP_400_BAD_REQUEST)
    
    @action(detail=True, methods=['get'])
    def status(self, request, pk=None):
        """获取计算状态"""
        try:
            job = Job.objects.get(pk=pk)
            return Response({
                'job_id': job.id,
                'status': job.status,
                'created_at': job.created_at,
                'updated_at': job.updated_at
            })
        except Job.DoesNotExist:
            return Response({'error': '作业不存在'}, status=status.HTTP_404_NOT_FOUND)
    
    @action(detail=True, methods=['get'])
    def results(self, request, pk=None):
        """获取计算结果"""
        try:
            job = Job.objects.get(pk=pk)
            if job.status != 'completed':
                return Response({
                    'job_id': job.id,
                    'status': job.status,
                    'message': '计算尚未完成'
                })
            
            try:
                result = job.result
                serializer = SimulationResultSerializer(result)
                return Response(serializer.data)
            except SimulationResult.DoesNotExist:
                return Response({'error': '结果不存在'}, status=status.HTTP_404_NOT_FOUND)
        except Job.DoesNotExist:
            return Response({'error': '作业不存在'}, status=status.HTTP_404_NOT_FOUND)
```

## 五、Celery任务设计

```python
# tasks.py
from celery import shared_task
from .models import Job, ElectrolyteFormula, SimulationResult
from molyte_cursor.src.io.file_generator import FileGenerator
from molyte_cursor.src.utils.command_executor import CommandExecutor
from molyte_cursor.src.analysis import Analyzer
import os
import json
import tempfile
import shutil

@shared_task
def run_electrolyte_simulation(job_id):
    """运行电解液模拟计算的Celery任务"""
    try:
        # 获取作业和配方数据
        job = Job.objects.get(pk=job_id)
        job.status = 'running'
        job.save()
        
        formula = ElectrolyteFormula.objects.get(job=job)
        
        # 创建临时工作目录
        work_dir = tempfile.mkdtemp(prefix=f"electrolyte_{job_id}_")
        
        try:
            # 1. 生成LAMMPS输入文件
            file_generator = FileGenerator(work_dir=work_dir)
            input_files = file_generator.generate_lammps_inputs(
                cations=formula.cations,
                anions=formula.anions,
                solvents=formula.solvents,
                temperatures=formula.temperatures,
                concentrations=formula.concentrations,
                calculation_types=formula.calculation_types,
                parameters=formula.simulation_parameters
            )
            
            # 2. 执行LAMMPS计算
            executor = CommandExecutor()
            success, log_file = executor.run_lammps(
                input_file=input_files['main_input'],
                work_dir=work_dir
            )
            
            if not success:
                raise Exception(f"LAMMPS执行失败，详见日志: {log_file}")
            
            # 3. 分析结果
            analyzer = Analyzer(
                base_dir=work_dir,
                output_dir=os.path.join(work_dir, 'analysis')
            )
            
            analysis_results = {}
            result_files = {}
            plots = {}
            
            # 根据计算类型执行不同分析
            if 'single_molecule' in formula.calculation_types:
                molecule_results = analyzer.run_gaussian_analysis(
                    directory=os.path.join(work_dir, 'gaussian'),
                    output_name='molecule_properties'
                )
                analysis_results['molecule'] = molecule_results
                
            if 'solvation' in formula.calculation_types:
                rdf_results = analyzer.run_rdf_analysis(job.name)
                analysis_results['rdf'] = rdf_results
                
            if 'transport' in formula.calculation_types:
                msd_results = analyzer.run_msd_analysis(job.name)
                electrochem_results = analyzer.run_electrochemical_analysis(
                    sample_name=job.name,
                    temperature=formula.temperatures[0]['value'],
                    concentration=formula.concentrations[0]['value'],
                    cation_valence=1,
                    anion_valence=1
                )
                analysis_results['msd'] = msd_results
                analysis_results['electrochemical'] = electrochem_results
            
            # 复制结果文件到永久存储位置
            # ...（这里需要实现文件复制逻辑）
            
            # 4. 保存结果
            result = SimulationResult.objects.create(
                job=job,
                result_files=result_files,
                analysis_results=analysis_results,
                plots=plots,
                summary=f"电解液{job.name}的模拟计算已完成"
            )
            
            job.status = 'completed'
            job.save()
            
            return f"作业 {job_id} 成功完成"
            
        except Exception as e:
            job.status = 'failed'
            job.save()
            raise e
        
        finally:
            # 清理临时文件
            # shutil.rmtree(work_dir)  # 在生产环境可以删除临时文件
            pass
            
    except Job.DoesNotExist:
        return f"作业 {job_id} 不存在"
    except ElectrolyteFormula.DoesNotExist:
        job.status = 'failed'
        job.save()
        return f"作业 {job_id} 的配方数据不存在"
    except Exception as e:
        try:
            job = Job.objects.get(pk=job_id)
            job.status = 'failed'
            job.save()
        except:
            pass
        return f"作业 {job_id} 执行出错: {str(e)}"
```

## 六、文件生成器设计

```python
# 与Molyte集成的LAMMPS文件生成器
class LammpsFileGenerator:
    """LAMMPS输入文件生成器"""
    
    def __init__(self, work_dir):
        self.work_dir = work_dir
        
    def generate_input_files(self, cations, anions, solvents, temperatures, 
                            concentrations, calculation_types, parameters=None):
        """生成LAMMPS输入文件"""
        # 实现基于电解液配方生成LAMMPS输入文件的逻辑
        # ...
        
        result_files = {
            'main_input': os.path.join(self.work_dir, 'electrolyte.in'),
            'data_file': os.path.join(self.work_dir, 'electrolyte.data'),
            'molecule_files': [...],
            'config_files': [...]
        }
        
        return result_files
```

## 七、Excel模板处理器

```python
# 处理Excel模板上传和下载的类
class ExcelTemplateHandler:
    """处理Excel模板上传和下载"""
    
    def __init__(self):
        pass
        
    def generate_template(self, output_path):
        """生成Excel模板文件"""
        # 使用pandas或openpyxl生成模板
        # ...
        
    def parse_template(self, input_file):
        """解析上传的Excel文件内容"""
        # 使用pandas读取Excel并解析为API格式数据
        # ...
        
        return {
            'cations': [...],
            'anions': [...],
            'solvents': [...],
            'temperatures': [...],
            'concentrations': [...]
        }
```

## 八、实现路线图

### 1. 项目初始化与基础架构设置（1-2周）
- 创建Django项目和应用
- 设置数据库和模型
- 配置Celery和Redis
- 实现基本API骨架
- 设置认证和权限

### 2. Molyte集成与文件生成功能（2-3周）
- 集成Molyte库
- 开发LAMMPS输入文件生成器
- 实现命令执行器和日志记录
- 实现Excel模板处理功能

### 3. 计算任务管理系统（2周）
- 实现Celery任务定义
- 开发作业状态管理
- 设计任务队列和监控机制
- 实现错误处理和通知系统

### 4. 结果分析与可视化（2-3周）
- 集成Molyte分析模块
- 实现RDF、MSD等分析功能
- 开发结果数据处理和存储
- 设计图表生成和导出功能

### 5. API完善与前端集成（1-2周）
- 完善所有API端点实现
- 编写API文档和测试
- 与前端进行集成测试
- 优化API性能和响应

### 6. 测试、优化与部署（2周）
- 编写单元测试和集成测试
- 进行性能测试和优化
- 准备部署配置和文档
- 部署到生产环境

## 九、技术挑战与解决方案

### 1. 长时间运行计算任务
**挑战**：分子动力学模拟可能需要长时间运行
**解决方案**：
- 使用Celery任务队列处理异步计算
- 实现计算状态跟踪和断点续传
- 为长时间任务设置超时和资源限制

### 2. 计算资源管理
**挑战**：多用户并发计算可能导致资源竞争
**解决方案**：
- 实现计算资源分配策略
- 设置任务优先级和队列管理
- 考虑使用容器化技术隔离计算环境

### 3. 结果数据存储与检索
**挑战**：计算结果可能非常大
**解决方案**：
- 设计高效的数据存储方案
- 使用文件系统或对象存储处理大文件
- 实现数据分页和懒加载机制

### 4. 系统安全性
**挑战**：运行用户提交的计算可能存在安全风险
**解决方案**：
- 严格验证和清理用户输入
- 使用沙盒环境执行计算
- 限制计算资源和访问权限

## 十、未来扩展方向

1. **批处理功能**：支持批量提交多组电解液配方进行对比研究
2. **自定义分析**：允许用户定义自己的分析方法和可视化
3. **机器学习集成**：添加机器学习预测功能，如电导率和稳定性预测
4. **实验数据比较**：支持上传实验数据与模拟结果进行比较
5. **高级可视化**：实现3D分子结构可视化和动态模拟轨迹查看 