#!/usr/bin/env python3
"""
命令行工具，用于处理从网页提交的电解液配方
"""
import os
import sys
import json
import argparse
from pathlib import Path
from io.recipe_processor import RecipeProcessor
from io.electrolyte_file_generator import generate_electrolyte_input_files
from utils.logger import Logger

def main():
    """主函数"""
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='处理电解液配方并生成输入文件')
    parser.add_argument('--input', '-i', required=True, help='输入的JSON配方文件路径')
    parser.add_argument('--output-dir', '-o', required=True, help='输出目录路径')
    args = parser.parse_args()
    
    # 初始化日志记录器
    logger = Logger().get_logger()
    logger.info("开始处理电解液配方")
    
    try:
        # 读取JSON配方文件
        with open(args.input, 'r', encoding='utf-8') as f:
            recipe_data = json.load(f)
        
        # 初始化配方处理器
        processor = RecipeProcessor()
        
        # 加载并验证配方
        recipe = processor.load_recipe_from_json(recipe_data)
        if not processor.validate_recipe(recipe):
            logger.error("配方验证失败")
            sys.exit(1)
        
        # 创建输出目录
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成所有必需的输入文件
        result = generate_electrolyte_input_files(recipe, str(output_dir))
        
        # 保存生成结果到JSON文件
        result_file = output_dir / 'generation_result.json'
        with open(result_file, 'w', encoding='utf-8') as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
        
        logger.info(f"电解液输入文件生成完成，结果保存在: {result_file}")
        
    except Exception as e:
        logger.error(f"处理配方时出错: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 