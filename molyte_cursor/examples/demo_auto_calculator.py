#!/usr/bin/env python3
"""
分子计算器演示脚本

该脚本展示了如何使用分子计算器自动计算模拟所需的分子数量，
而不是依赖于Excel中手动计算的结果。
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

def demo_molecule_calculator():
    """展示分子计算器的基本用法"""
    logger = setup_logger('calculator_demo', console_level=logging.INFO)
    
    print("\n===== 分子计算器基本用法演示 =====")
    
    # 初始化分子计算器
    calculator = MoleculeCalculator(logger=logger)
    
    # 计算简单的离子和溶剂数量
    print("\n1. 简单计算 - 给定体积、浓度和溶剂比例")
    result = calculator.calculate_molecule_counts(
        volume=50*50*50,           # 体积（立方埃）
        concentration=1.0,         # 浓度（mol/L）
        solvent_ratio=10.0,        # 溶剂比例（相对于离子对）
        volume_unit='angstrom3',   # 体积单位
        concentration_unit='mol/L' # 浓度单位
    )
    
    print(f"计算结果:")
    print(f"  阳离子数量: {result['cation']}")
    print(f"  阴离子数量: {result['anion']}")
    print(f"  溶剂分子数量: {result['solvent']}")
    print(f"  总分子数量: {result['total']}")
    
    # 从系统定义计算
    print("\n2. 从系统定义计算 - 更复杂的系统")
    system_result = calculator.calculate_from_system_definition(
        box_size=(50, 50, 50),         # 盒子尺寸（埃）
        cation_concentration=1.0,      # 阳离子浓度（mol/L）
        solvent_data={                 # 溶剂数据（名称：比例）
            'solvent1': 7.0,           # 70%的溶剂1
            'solvent2': 3.0            # 30%的溶剂2
        },
        cation_types=2,                # 2种阳离子
        anion_types=1,                 # 1种阴离子
        length_unit='angstrom',        # 长度单位
        concentration_unit='mol/L'     # 浓度单位
    )
    
    print(f"系统计算结果:")
    print(f"  阳离子: {system_result['cations']}")
    print(f"  阴离子: {system_result['anions']}")
    print(f"  溶剂: {system_result['solvents']}")
    print(f"  总分子数量: {system_result['total']}")
    
    # 自定义组分计算
    print("\n3. 自定义组分计算")
    custom_result = calculator.calculate_custom_composition(
        box_size=(50, 50, 50),
        components={
            'Li': {
                'type': 'cation',
                'concentration': 0.8
            },
            'Na': {
                'type': 'cation',
                'concentration': 0.2
            },
            'TFSI': {
                'type': 'anion',
                'concentration': 1.0
            },
            'EC': {
                'type': 'solvent',
                'ratio': 7.0  # 相对于离子对的比例
            },
            'DMC': {
                'type': 'solvent',
                'ratio': 3.0  # 相对于离子对的比例
            }
        }
    )
    
    print(f"自定义组分计算结果:")
    print(f"  阳离子: {custom_result['cations']}")
    print(f"  阴离子: {custom_result['anions']}")
    print(f"  溶剂: {custom_result['solvents']}")
    print(f"  总分子数量: {custom_result['total']}")

def create_demo_excel():
    """创建演示用的Excel模板"""
    output_file = os.path.join(os.path.dirname(__file__), 'templates', 'auto_calculator_template.xlsx')
    
    # 创建DataFrame
    data = {
        'name': ['System1', 'System2', 'System3', 'System4'],
        
        # 盒子尺寸
        'box_size_x': [50, 60, 50, 70],
        'box_size_y': [50, 60, 50, 70],
        'box_size_z': [50, 60, 50, 70],
        
        # 浓度和比例
        'concentration': [1.0, 0.5, 1.0, 2.0],
        'solvent_ratio': [10.0, 20.0, 15.0, 5.0],
        'cation_anion_ratio': [1.0, 1.0, 2.0, 0.5],
        
        # 阳离子、阴离子和溶剂列（空值，将由自动计算填充）
        'cation1': [None, None, None, 150],  # System4手动指定数量
        'cation2': [None, None, None, None],
        'anion1': [None, None, None, 300],   # System4手动指定数量
        'sol1': [None, None, None, None],
        'sol2': [None, None, None, None],
        
        # 分子名称和SMILE
        'cation1_name': ['Li', 'Na', 'Li', 'K'],
        'cation2_name': ['Na', 'K', 'Na', 'Li'],
        'anion1_name': ['TFSI', 'PF6', 'BF4', 'Cl'],
        'sol1_name': ['EC', 'DMC', 'EC', 'DMSO'],
        'sol2_name': ['DMC', 'PC', 'DMC', 'H2O']
    }
    
    df = pd.DataFrame(data)
    
    # 保存到Excel
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    df.to_excel(output_file, index=False)
    
    print(f"\n已创建演示Excel模板: {output_file}")
    return output_file

def demo_excel_processor(excel_file):
    """展示Excel处理器的使用"""
    logger = setup_logger('excel_processor_demo', console_level=logging.INFO)
    
    print("\n===== Excel处理器演示 =====")
    
    # 初始化Excel处理器
    processor = ExcelProcessor(logger=logger)
    
    # 处理Excel文件
    processed_file = os.path.join(os.path.dirname(excel_file), 'auto_calculator_processed.xlsx')
    result_df = processor.process_excel_file(excel_file, output_file=processed_file)
    
    print(f"\n已处理Excel文件并保存到: {processed_file}")
    
    # 打印处理结果
    print("\n处理后的分子数量:")
    for idx, row in result_df.iterrows():
        print(f"\n系统: {row['name']}")
        print(f"  阳离子: cation1={row['cation1']}, cation2={row['cation2']}")
        print(f"  阴离子: anion1={row['anion1']}")
        print(f"  溶剂: sol1={row['sol1']}, sol2={row['sol2']}")

def demo_excel_reader(excel_file):
    """展示带有自动计算功能的Excel读取器的使用"""
    logger = setup_logger('excel_reader_demo', console_level=logging.INFO)
    
    print("\n===== Excel读取器演示（带自动计算） =====")
    
    # 初始化Excel读取器，启用自动计算
    reader = ExcelReader(auto_calculate=True)
    
    # 读取Excel文件
    systems = reader.read_excel(excel_file)
    
    # 打印系统信息
    print(f"\n共读取 {len(systems)} 个分子系统:")
    for system in systems:
        print(f"\n系统: {system.name}")
        print(f"  阳离子: {[f'{c.name}({c.count})' for c in system.cations]}")
        print(f"  阴离子: {[f'{a.name}({a.count})' for a in system.anions]}")
        print(f"  溶剂: {[f'{s.name}({s.count})' for s in system.solvents]}")

if __name__ == "__main__":
    # 演示分子计算器
    demo_molecule_calculator()
    
    # 创建演示Excel文件
    excel_file = create_demo_excel()
    
    # 演示Excel处理器
    demo_excel_processor(excel_file)
    
    # 演示Excel读取器
    demo_excel_reader(excel_file) 