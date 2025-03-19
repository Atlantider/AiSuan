import sys
import os
from typing import Dict, Any

# 打印Python路径
print("Python路径:")
for path in sys.path:
    print(f"  - {path}")

# 添加molyte_cursor到路径
MOLYTE_PATH = os.path.join(os.getcwd(), 'molyte_cursor')
if MOLYTE_PATH not in sys.path:
    sys.path.append(MOLYTE_PATH)
    print(f"\n已添加 {MOLYTE_PATH} 到Python路径")

# 尝试导入模块
print("\n尝试导入molyte_cursor模块...")
try:
    from molyte_cursor.src.core.main import generate_input_files
    from molyte_cursor.src.io.file_generator import LAMMPSFileGenerator
    print("成功导入molyte_cursor模块!")
    
    # 测试不同配置
    test_configs = [
        {
            "name": "标准配方",
            "config": {
                "formulation_name": "标准LiPF6/EC配方",
                "temperature": 300,
                "box_size": 40,
                "concentration": 1.0,
                "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
                "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
            },
            "output_dir": os.path.join(os.getcwd(), "test_output/standard")
        },
        {
            "name": "高温配方",
            "config": {
                "formulation_name": "高温LiPF6/EC配方",
                "temperature": 350,
                "box_size": 40,
                "concentration": 1.0,
                "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
                "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
            },
            "output_dir": os.path.join(os.getcwd(), "test_output/high_temp")
        },
        {
            "name": "混合溶剂配方",
            "config": {
                "formulation_name": "LiPF6/EC-DMC混合溶剂配方",
                "temperature": 300,
                "box_size": 40,
                "concentration": 1.0,
                "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
                "solvents": [
                    {"name": "EC", "smile": "C1OC(=O)O1", "concentration": 0.5},
                    {"name": "DMC", "smile": "COC(=O)OC", "concentration": 0.5}
                ]
            },
            "output_dir": os.path.join(os.getcwd(), "test_output/mixed_solvent")
        },
        {
            "name": "混合盐配方",
            "config": {
                "formulation_name": "LiPF6-LiBF4混合盐配方",
                "temperature": 300,
                "box_size": 40,
                "concentration": 1.0,
                "salts": [
                    {"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 0.7},
                    {"name": "LiBF4", "cation": "Li", "anion": "BF4", "concentration": 0.3}
                ],
                "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
            },
            "output_dir": os.path.join(os.getcwd(), "test_output/mixed_salt")
        }
    ]
    
    generator = LAMMPSFileGenerator()
    print("\n创建LAMMPSFileGenerator实例成功")
    
    # 测试每个配置
    for test in test_configs:
        print(f"\n测试配置: {test['name']}")
        print(f"输出目录: {test['output_dir']}")
        
        # 使用generate_input_files函数
        result = generate_input_files(test['config'], test['output_dir'])
        print(f"生成文件结果: {result}")
        
        # 验证文件是否存在
        inp_file = result['inp_file']
        lammps_file = result['lammps_file']
        
        print(f"INP文件存在: {os.path.exists(inp_file)}")
        print(f"LAMMPS文件存在: {os.path.exists(lammps_file)}")
    
except ImportError as e:
    print(f"导入错误: {e}")
except Exception as e:
    print(f"其他错误: {e}") 