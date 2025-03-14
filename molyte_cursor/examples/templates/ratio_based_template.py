#!/usr/bin/env python3
"""
基于比例计算的Excel模板生成器

该脚本生成一个Excel模板，展示如何使用基于主阳离子和比例的计算方式
"""

import os
import sys
import pandas as pd

# 添加项目根目录到路径
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
sys.path.insert(0, project_root)

def create_ratio_based_template():
    """创建基于比例的Excel模板"""
    
    # 定义输出文件路径
    output_file = os.path.join(script_dir, 'ratio_based_template.xlsx')
    
    # 创建示例数据
    data = {
        # 基本信息
        'name': ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
        
        # 立方体盒子尺寸
        'box_size': [50, 60, 70, 80, 90],
        
        # 主阳离子浓度
        'concentration': [1.0, 0.5, 1.5, 2.0, 0.8],
        
        # 阳离子列（将由计算填充）
        'cation1': [None, None, None, None, None],
        'cation2': [None, None, None, None, None],
        
        # 阴离子列（将由计算填充）
        'anion1': [None, None, None, None, None],
        'anion2': [None, None, None, None, None],
        
        # 溶剂列（将由计算填充）
        'sol1': [None, None, None, None, None],
        'sol2': [None, None, None, None, None],
        'sol3': [None, None, None, None, None],
        
        # 比例设置 - 以cation1为基准(1.0)
        'cation1_ratio': [1.0, 1.0, 1.0, 1.0, 1.0],  # 主阳离子，比例始终为1.0
        'cation2_ratio': [0.5, 0.3, 0.2, 0.1, 0.4],  # 第二种阳离子的比例
        'anion1_ratio': [1.5, 1.0, 1.2, 2.0, 1.3],   # 第一种阴离子的比例
        'anion2_ratio': [0.0, 0.3, 0.0, 0.1, 0.1],   # 第二种阴离子的比例
        'sol1_ratio': [5.0, 10.0, 3.0, 8.0, 6.0],    # 第一种溶剂的比例
        'sol2_ratio': [5.0, 0.0, 3.0, 0.0, 4.0],     # 第二种溶剂的比例
        'sol3_ratio': [0.0, 0.0, 3.0, 2.0, 0.0],     # 第三种溶剂的比例
        
        # 分子名称
        'cation1_name': ['Li', 'Na', 'Li', 'K', 'Li'],
        'cation2_name': ['Na', 'K', 'Na', 'Li', 'Mg'],
        'anion1_name': ['TFSI', 'PF6', 'BF4', 'Cl', 'FSI'],
        'anion2_name': ['PF6', 'Cl', 'I', 'Br', 'TFSI'],
        'sol1_name': ['EC', 'DMC', 'PC', 'DMSO', 'EC'],
        'sol2_name': ['DMC', 'EC', 'DMC', 'PC', 'DMC'],
        'sol3_name': ['PC', 'PC', 'FEC', 'H2O', 'DEC'],
        
        # SMILE字符串（可选）
        'cation1_smile': ['[Li+]', '[Na+]', '[Li+]', '[K+]', '[Li+]'],
        'anion1_smile': ['N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F', '[PF6-]', '[BF4-]', '[Cl-]', '[N-]S(=O)(=O)F'],
        
        # 其他参数（可选）
        'temperature': [300, 298, 310, 350, 320],
        'density': [None, None, None, None, None]  # 可以通过计算密度来验证
    }
    
    # 创建DataFrame
    df = pd.DataFrame(data)
    
    # 保存到Excel
    df.to_excel(output_file, index=False)
    
    print(f"已创建基于比例的Excel模板: {output_file}")
    return output_file

if __name__ == "__main__":
    create_ratio_based_template() 