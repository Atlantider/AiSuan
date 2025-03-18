#!/usr/bin/env python
"""
测试脚本 - 在Django环境中测试molyte_cursor模块的功能
"""
import os
import sys
import django

# 添加项目根目录到Python路径
sys.path.append('/Users/xiaoji/Documents/AiSuan')
sys.path.append('/Users/xiaoji/Documents/AiSuan/backend')

# 设置Django环境
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'electrolyte_backend.settings')
django.setup()

print("Django环境已设置完成")

# 尝试导入molyte_cursor模块
try:
    from molyte_cursor.src.core.main import generate_input_files
    from molyte_cursor.src.io.file_generator import LAMMPSFileGenerator
    print("成功导入molyte_cursor模块!")
    
    # 测试配置
    config = {
        "formulation_name": "在Django中测试的配方",
        "temperature": 300,
        "box_size": 40,
        "concentration": 1.0,
        "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
        "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
    }
    
    # 使用Django项目中的media目录作为输出
    from django.conf import settings
    output_dir = os.path.join(settings.MEDIA_ROOT, 'test_output')
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"输出目录: {output_dir}")
    
    # 创建LAMMPSFileGenerator实例
    generator = LAMMPSFileGenerator()
    print("创建LAMMPSFileGenerator实例成功")
    
    # 生成文件
    result = generate_input_files(config, output_dir)
    print(f"生成文件结果: {result}")
    print(f"INP文件路径: {result['inp_file']}")
    print(f"LAMMPS文件路径: {result['lammps_file']}")
    
    # 验证文件是否存在
    print(f"INP文件存在: {os.path.exists(result['inp_file'])}")
    print(f"LAMMPS文件存在: {os.path.exists(result['lammps_file'])}")
    
    # 读取并打印文件内容
    print("\nINP文件内容:")
    with open(result['inp_file'], 'r') as f:
        print(f.read())
    
    print("\nLAMMPS文件内容:")
    with open(result['lammps_file'], 'r') as f:
        print(f.read())
        
    print("\n测试完成！")
    
except ImportError as e:
    print(f"导入错误: {e}")
except Exception as e:
    print(f"其他错误: {e}") 