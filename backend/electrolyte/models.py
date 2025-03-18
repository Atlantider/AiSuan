from django.db import models
from django.contrib.auth.models import User
import uuid
import os

def calculation_file_path(instance, filename):
    """为计算文件生成唯一的文件路径"""
    ext = filename.split('.')[-1]
    filename = f"{uuid.uuid4()}.{ext}"
    return os.path.join('calculations', str(instance.id), filename)

class Solvent(models.Model):
    """溶剂模型"""
    name = models.CharField(max_length=100, verbose_name="溶剂名称")
    smile = models.CharField(max_length=255, verbose_name="SMILE表示法")
    description = models.TextField(blank=True, null=True, verbose_name="描述")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "溶剂"
        verbose_name_plural = "溶剂"

class Salt(models.Model):
    """盐类模型"""
    name = models.CharField(max_length=100, verbose_name="盐名称")
    cation = models.CharField(max_length=100, verbose_name="阳离子")
    anion = models.CharField(max_length=100, verbose_name="阴离子")
    description = models.TextField(blank=True, null=True, verbose_name="描述")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "盐"
        verbose_name_plural = "盐"

class ElectrolyteFormulation(models.Model):
    """电解质配方模型"""
    name = models.CharField(max_length=255, verbose_name="配方名称")
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="formulations", verbose_name="用户")
    description = models.TextField(blank=True, null=True, verbose_name="描述")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = "电解质配方"
        verbose_name_plural = "电解质配方"

class FormulationComponent(models.Model):
    """配方组件模型（盐和溶剂）"""
    COMPONENT_TYPES = [
        ('salt', '盐'),
        ('solvent', '溶剂'),
    ]
    
    formulation = models.ForeignKey(ElectrolyteFormulation, on_delete=models.CASCADE, related_name="components", verbose_name="配方")
    component_type = models.CharField(max_length=10, choices=COMPONENT_TYPES, verbose_name="组件类型")
    salt = models.ForeignKey(Salt, on_delete=models.CASCADE, null=True, blank=True, verbose_name="盐")
    solvent = models.ForeignKey(Solvent, on_delete=models.CASCADE, null=True, blank=True, verbose_name="溶剂")
    concentration = models.FloatField(verbose_name="浓度/比例")
    
    def __str__(self):
        if self.component_type == 'salt':
            return f"{self.salt.name} ({self.concentration} mol/L)"
        else:
            return f"{self.solvent.name} ({self.concentration}%)"
    
    class Meta:
        verbose_name = "配方组件"
        verbose_name_plural = "配方组件"

class SimulationParameters(models.Model):
    """模拟参数模型"""
    formulation = models.OneToOneField(ElectrolyteFormulation, on_delete=models.CASCADE, related_name="parameters", verbose_name="配方")
    temperature = models.FloatField(verbose_name="温度(K)")
    pressure = models.FloatField(verbose_name="压力(atm)", default=1.0)
    time_step = models.FloatField(verbose_name="时间步长(fs)", default=1.0)
    equilibration_steps = models.IntegerField(verbose_name="平衡步数", default=10000)
    production_steps = models.IntegerField(verbose_name="生产步数", default=50000)
    cutoff = models.FloatField(verbose_name="截断距离(Å)", default=12.0)
    additional_params = models.JSONField(blank=True, null=True, verbose_name="额外参数")
    
    def __str__(self):
        return f"{self.formulation.name} - {self.temperature}K"
    
    class Meta:
        verbose_name = "模拟参数"
        verbose_name_plural = "模拟参数"

class Calculation(models.Model):
    """计算模型"""
    STATUS_CHOICES = [
        ('pending', '待处理'),
        ('running', '运行中'),
        ('completed', '已完成'),
        ('failed', '失败'),
    ]
    
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name="calculations", verbose_name="用户")
    formulation = models.ForeignKey(ElectrolyteFormulation, on_delete=models.CASCADE, related_name="calculations", verbose_name="配方")
    name = models.CharField(max_length=255, verbose_name="计算名称")
    description = models.TextField(blank=True, null=True, verbose_name="描述")
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending', verbose_name="状态")
    task_id = models.CharField(max_length=255, blank=True, null=True, verbose_name="Celery任务ID")
    input_file = models.FileField(upload_to=calculation_file_path, blank=True, null=True, verbose_name="输入文件")
    output_file = models.FileField(upload_to=calculation_file_path, blank=True, null=True, verbose_name="输出文件")
    log_file = models.FileField(upload_to=calculation_file_path, blank=True, null=True, verbose_name="日志文件")
    error_message = models.TextField(blank=True, null=True, verbose_name="错误信息")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")
    started_at = models.DateTimeField(blank=True, null=True, verbose_name="开始时间")
    finished_at = models.DateTimeField(blank=True, null=True, verbose_name="完成时间")
    
    def __str__(self):
        return self.name
    
    class Meta:
        verbose_name = "计算"
        verbose_name_plural = "计算"

class CalculationResult(models.Model):
    """计算结果模型"""
    calculation = models.OneToOneField(Calculation, on_delete=models.CASCADE, related_name="result", verbose_name="计算")
    ionic_conductivity = models.FloatField(blank=True, null=True, verbose_name="离子电导率(S/m)")
    diffusion_coefficients = models.JSONField(blank=True, null=True, verbose_name="扩散系数")
    radial_distribution = models.JSONField(blank=True, null=True, verbose_name="径向分布函数")
    density = models.FloatField(blank=True, null=True, verbose_name="密度(g/cm³)")
    viscosity = models.FloatField(blank=True, null=True, verbose_name="黏度(cP)")
    additional_results = models.JSONField(blank=True, null=True, verbose_name="额外结果")
    
    def __str__(self):
        return f"{self.calculation.name} - 结果"
    
    class Meta:
        verbose_name = "计算结果"
        verbose_name_plural = "计算结果"

class InputFile(models.Model):
    """INP输入文件模型"""
    formulation = models.OneToOneField(ElectrolyteFormulation, on_delete=models.CASCADE, related_name="input_file", verbose_name="配方")
    content = models.TextField(verbose_name="文件内容")
    file_path = models.CharField(max_length=255, blank=True, null=True, verbose_name="文件路径")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")
    
    def save(self, *args, **kwargs):
        # 保存到数据库时，同时写入文件系统
        super().save(*args, **kwargs)
        
        # 确保目录存在
        from django.conf import settings
        import os
        
        file_dir = os.path.join(settings.MEDIA_ROOT, 'input_files')
        os.makedirs(file_dir, exist_ok=True)
        
        # 生成文件路径
        file_path = os.path.join(file_dir, f'electrolyte_input_{self.formulation.id}.inp')
        
        # 写入文件
        with open(file_path, 'w') as f:
            f.write(self.content)
        
        # 更新文件路径并再次保存
        if self.file_path != file_path:
            self.file_path = file_path
            super().save(update_fields=['file_path'])
    
    def __str__(self):
        return f"输入文件 - {self.formulation.name}"
    
    class Meta:
        verbose_name = "输入文件"
        verbose_name_plural = "输入文件"
