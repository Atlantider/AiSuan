"""
电解液模块辅助功能
"""

import os
import logging
import json
from django.conf import settings
from django.db import models
from .models import ElectrolyteFormulation

logger = logging.getLogger(__name__)

# 添加一个函数直接修复指定配方的SMILE字符串
def fix_smile_in_formulation(formulation_id):
    """
    修复特定配方中的SMILE字符串，尤其是EC的错误SMILE字符串。
    
    Args:
        formulation_id: 配方的ID
        
    Returns:
        tuple: (成功标志, 消息)
    """
    try:
        # 获取配方对象
        formulation = ElectrolyteFormulation.objects.get(id=formulation_id)
    except ElectrolyteFormulation.DoesNotExist:
        logger.error(f"找不到ID为{formulation_id}的配方")
        return False, f"找不到ID为{formulation_id}的配方"
    
    # 获取输入文件路径
    input_file_path = None
    if formulation.input_file:
        input_file_path = formulation.input_file.path
    
    if not input_file_path or not os.path.exists(input_file_path):
        logger.warning(f"配方ID {formulation_id} 没有输入文件或文件不存在")
        return False, f"配方ID {formulation_id} 没有输入文件或文件不存在"
    
    # 读取输入文件内容
    modified = False
    try:
        with open(input_file_path, 'r') as f:
            content = f.read()
            
        # 检查并修复EC的SMILE字符串
        if "O)O1" in content:
            logger.info(f"在配方ID {formulation_id} 中发现错误的SMILE字符串: O)O1")
            content = content.replace("O)O1", "C1COC(=O)O1")
            
            # 写回修复后的内容
            with open(input_file_path, 'w') as f:
                f.write(content)
                
            logger.info(f"成功修复配方ID {formulation_id} 中的SMILE字符串")
            modified = True
        else:
            logger.info(f"配方ID {formulation_id} 中没有发现错误的SMILE字符串")
    
    except Exception as e:
        error_msg = f"修复配方ID {formulation_id} 的SMILE字符串时出错: {str(e)}"
        logger.error(error_msg)
        return False, error_msg
    
    if modified:
        return True, f"成功修复配方ID {formulation_id} 中的SMILE字符串"
    else:
        return True, f"配方ID {formulation_id} 中没有需要修复的SMILE字符串"

# 批量修复所有EC的错误SMILE字符串
def fix_all_ec_smile_errors():
    """
    修复所有配方中EC的错误SMILE字符串。
    
    Returns:
        tuple: (修复成功的数量, 修复失败的数量, 详细信息列表)
    """
    formulations = ElectrolyteFormulation.objects.all()
    fixed_count = 0
    error_count = 0
    details = []
    
    for formulation in formulations:
        success, message = fix_smile_in_formulation(formulation.id)
        
        details.append({
            "formulation_id": formulation.id,
            "success": success,
            "message": message
        })
        
        if success and "成功修复" in message:
            fixed_count += 1
        elif not success:
            error_count += 1
    
    logger.info(f"SMILE字符串修复任务完成：成功修复 {fixed_count} 个配方，{error_count} 个错误")
    return fixed_count, error_count, details 