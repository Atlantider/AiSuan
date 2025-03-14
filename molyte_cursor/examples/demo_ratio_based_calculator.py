#!/usr/bin/env python3
"""
基于比例的分子计算器演示脚本

该脚本演示了如何使用基于主阳离子和比例的分子计算器，
相比原来的自动计算方式，这种计算方式能更精确地控制各种组分的比例。
"""

import os
import sys
import logging
import pandas as pd

# 添加项目根目录到路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from molyte_cursor.src.utils.molecule_calculator import MoleculeCalculator
from molyte_cursor.src.io.excel_processor import ExcelProcessor
from molyte_cursor.src.io.excel_reader import ExcelReader
from molyte_cursor.src.utils.logger import setup_logger
from molyte_cursor.examples.templates.ratio_based_template import create_ratio_based_template

def demo_ratio_based_calculator():
    """展示基于比例的分子计算器的基本用法"""
    logger = setup_logger('ratio_calculator_demo', console_level=logging.INFO)
    
    print("\n===== 基于主阳离子的比例计算演示 =====")
    
    # 初始化分子计算器
    calculator = MoleculeCalculator(logger=logger)
    
    # 示例 - 80×80×80 的立方体，主阳离子Li浓度1.0 mol/L
    # 其他组分按照如下比例计算：
    # Na（第二种阳离子）: Li = 0.5 : 1
    # TFSI（第一种阴离子）: Li = 1.5 : 1
    # PF6（第二种阴离子）: Li = 0 : 1 (不存在)
    # EC（第一种溶剂）: Li = 5 : 1
    # DMC（第二种溶剂）: Li = 5 : 1
    # PC（第三种溶剂）: Li = 0 : 1 (不存在)
    
    result = calculator.calculate_from_primary_cation(
        box_size=80,                      # 立方体盒子边长（埃）
        primary_cation_concentration=1.0, # 主阳离子浓度（mol/L）
        component_ratios={                # 各组分相对于主阳离子的比例
            'cation1': 1.0,               # Li (主阳离子，比例始终为1.0)
            'cation2': 0.5,               # Na
            'anion1': 1.5,                # TFSI
            'anion2': 0.0,                # PF6 (比例为0，不存在)
            'sol1': 5.0,                  # EC
            'sol2': 5.0,                  # DMC
            'sol3': 0.0                   # PC (比例为0，不存在)
        },
        length_unit='angstrom',
        concentration_unit='mol/L'
    )
    
    print("\n基于主阳离子的计算结果:")
    print(f"  阳离子: {result['cations']}")
    print(f"  阴离子: {result['anions']}")
    print(f"  溶剂: {result['solvents']}")
    print(f"  总分子数量: {result['total']}")
    
    return result

def generate_and_process_excel():
    """生成并处理基于比例的Excel模板"""
    logger = setup_logger('excel_processor_demo', console_level=logging.INFO)
    
    print("\n===== 基于比例的Excel处理演示 =====")
    
    # 创建模板Excel文件
    excel_file = create_ratio_based_template()
    
    # 初始化Excel处理器
    processor = ExcelProcessor(logger=logger)
    
    # 处理Excel文件
    processed_file = os.path.join(os.path.dirname(excel_file), 'ratio_based_processed.xlsx')
    result_df = processor.process_excel_file(excel_file, output_file=processed_file)
    
    print(f"\n已处理基于比例的Excel文件并保存到: {processed_file}")
    
    # 打印处理结果 - 只显示第一行作为示例
    row = result_df.iloc[0]
    print(f"\n样本 {row['name']} 的计算结果:")
    print(f"  盒子大小: {row['box_size']} 埃")
    print(f"  主阳离子浓度: {row['concentration']} mol/L")
    print(f"  阳离子: cation1({row['cation1']})={row['cation1_name']}, cation2({row['cation2']})={row['cation2_name']}")
    print(f"  阴离子: anion1({row['anion1']})={row['anion1_name']}, anion2({row['anion2']})={row['anion2_name']}")
    print(f"  溶剂: sol1({row['sol1']})={row['sol1_name']}, sol2({row['sol2']})={row['sol2_name']}, sol3({row['sol3']})={row['sol3_name']}")
    
    # 计算总分子数
    total = row['cation1'] + row['cation2'] + row['anion1'] + row['anion2'] + row['sol1'] + row['sol2'] + row['sol3']
    print(f"  总分子数量: {total}")
    
    # 打印比例
    print(f"\n组分比例:")
    print(f"  {row['cation1_name']}:{row['cation2_name']}:{row['anion1_name']}:{row['sol1_name']}:{row['sol2_name']} = 1.0:{row['cation2_ratio']}:{row['anion1_ratio']}:{row['sol1_ratio']}:{row['sol2_ratio']}")
    print(f"  实际计算结果 = 1.0:{row['cation2']/row['cation1']:.2f}:{row['anion1']/row['cation1']:.2f}:{row['sol1']/row['cation1']:.2f}:{row['sol2']/row['cation1']:.2f}")
    
    return result_df, excel_file

def demo_excel_reader(excel_file):
    """展示带有比例计算功能的Excel读取器的使用"""
    logger = setup_logger('excel_reader_demo', console_level=logging.INFO)
    
    print("\n===== 读取基于比例计算的Excel文件 =====")
    
    # 初始化Excel读取器，启用自动计算
    reader = ExcelReader(auto_calculate=True)
    
    # 读取Excel文件
    systems = reader.read_excel(excel_file)
    
    # 打印系统信息 - 只显示第一个系统作为示例
    system = systems[0]
    print(f"\n读取的系统: {system.name}")
    print(f"  阳离子: {[f'{c.name}({c.count})' for c in system.cations]}")
    print(f"  阴离子: {[f'{a.name}({a.count})' for a in system.anions]}")
    print(f"  溶剂: {[f'{s.name}({s.count})' for s in system.solvents]}")
    
    return systems

def calculate_example():
    """计算一个完整的示例"""
    # 设置模拟参数
    box_size = 70.0  # 立方体边长（埃）
    li_concentration = 1.0  # 主阳离子(Li)浓度（mol/L）
    
    # 定义组分比例 - 以Li为基准(1.0)
    component_ratios = {
        'cation1': 1.0,    # Li (主阳离子)
        'cation2': 0.0,    # Na (此例中不使用)
        'anion1': 1.0,     # TFSI 
        'sol1': 7.0,       # EC
        'sol2': 3.0        # DMC
    }
    
    # 计算分子数量
    calculator = MoleculeCalculator()
    result = calculator.calculate_from_primary_cation(
        box_size=box_size,
        primary_cation_concentration=li_concentration,
        component_ratios=component_ratios
    )
    
    # 打印结果
    print("\n===== 实际应用示例 =====")
    print(f"盒子: {box_size}×{box_size}×{box_size} 埃立方体")
    print(f"Li浓度: {li_concentration} mol/L")
    print(f"组分比例: Li:TFSI:EC:DMC = 1:1:7:3")
    print("\n计算结果:")
    print(f"  Li: {result['cations'].get('cation1', 0)} 个")
    print(f"  TFSI: {result['anions'].get('anion1', 0)} 个")
    print(f"  EC: {result['solvents'].get('sol1', 0)} 个")
    print(f"  DMC: {result['solvents'].get('sol2', 0)} 个")
    print(f"  总分子数: {result['total']} 个")
    
    return result

if __name__ == "__main__":
    # 展示基于比例的计算器
    demo_ratio_based_calculator()
    
    # 生成并处理基于比例的Excel
    result_df, excel_file = generate_and_process_excel()
    
    # 展示Excel读取器
    demo_excel_reader(excel_file)
    
    # 计算实际应用示例
    calculate_example() 