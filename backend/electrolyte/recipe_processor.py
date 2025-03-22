"""
电解液配方处理器，用于从网页接收配方并调用molyte_cursor处理
"""
import os
import sys
import json
import logging
import datetime
from pathlib import Path
from django.conf import settings
from .models import ElectrolyteFormulation, InputFile, Calculation, CalculationResult

# 配置日志
logger = logging.getLogger(__name__)

# 添加molyte_cursor到Python路径
MOLYTE_PATH = os.path.join(os.path.dirname(settings.BASE_DIR), 'molyte_cursor')
if MOLYTE_PATH not in sys.path:
    sys.path.append(MOLYTE_PATH)

# 尝试导入molyte_cursor功能
try:
    from molyte_cursor.src.io.recipe_processor import RecipeProcessor
    from molyte_cursor.src.io.electrolyte_file_generator import generate_electrolyte_input_files
    MOLYTE_AVAILABLE = True
    logger.info("成功导入molyte_cursor模块!")
except ImportError as e:
    logger.error(f"导入molyte_cursor错误: {e}")
    MOLYTE_AVAILABLE = False

def process_recipe_from_web(recipe_data, output_dir=None, user=None):
    """处理从网页提交的配方数据
    
    Args:
        recipe_data: 从前端提交的JSON格式配方数据
        output_dir: 指定输出目录，如果为None则自动创建
        user: 关联的用户对象
        
    Returns:
        字典，包含处理结果信息
    """
    logger.info(f"接收到配方数据: {recipe_data}")
    
    # 如果molyte_cursor模块不可用，返回错误
    if not MOLYTE_AVAILABLE:
        logger.error("molyte_cursor模块不可用，无法处理配方")
        return {
            "success": False,
            "error": "molyte_cursor模块不可用，请联系管理员"
        }
    
    try:
        # 初始化配方处理器
        processor = RecipeProcessor()
        
        # 加载并验证配方
        recipe = processor.load_recipe_from_json(recipe_data)
        if not processor.validate_recipe(recipe):
            logger.error("配方验证失败")
            return {
                "success": False,
                "error": "配方验证失败，请检查配方数据"
            }
        
        # 创建或获取电解液配方记录
        formulation, created = ElectrolyteFormulation.objects.get_or_create(
            name=recipe.name,
            defaults={
                "description": f"从网页提交的配方，温度: {recipe.temperature}K，盒子大小: {recipe.box_size}Å",
                "user": user,
                "created_at": datetime.datetime.now()
            }
        )
        
        if created:
            logger.info(f"创建了新的配方记录: {formulation.id}")
        else:
            logger.info(f"使用现有配方记录: {formulation.id}")
        
        # 创建输出目录
        if output_dir is None:
            output_dir = os.path.join(settings.MEDIA_ROOT, 'formulations', str(formulation.id))
        
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"输出目录: {output_dir}")
        
        # 生成所有必需的输入文件
        logger.info("开始生成输入文件...")
        result = generate_electrolyte_input_files(recipe, output_dir)
        logger.info(f"输入文件生成完成: {result}")
        
        # 创建或更新InputFile记录
        inp_file_path = os.path.join(output_dir, f"{recipe.name}.inp")
        
        # 创建或更新InputFile记录
        if hasattr(formulation, 'inputfile') and formulation.inputfile:
            input_file = formulation.inputfile
            input_file.file_path = inp_file_path
            input_file.save()
            logger.info(f"更新了现有InputFile记录: {input_file.id}")
        else:
            input_file = InputFile.objects.create(
                formulation=formulation,
                file_path=inp_file_path
            )
            logger.info(f"创建了新的InputFile记录: {input_file.id}")
        
        # 创建计算任务
        calculation = Calculation.objects.create(
            formulation=formulation,
            user=user,
            name=f"{recipe.name} 计算任务",
            description=f"从网页提交的配方 {recipe.name} 的计算任务",
            status='pending'
        )
        logger.info(f"创建了新的计算任务: {calculation.id}")
        
        # 返回成功结果
        return {
            "success": True,
            "formulation_id": formulation.id,
            "calculation_id": calculation.id,
            "output_dir": output_dir,
            "generated_files": result
        }
        
    except Exception as e:
        logger.error(f"处理配方时出错: {e}", exc_info=True)
        return {
            "success": False,
            "error": f"处理配方时出错: {str(e)}"
        } 