import os
import sys
import django
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='处理电解液配方INP文件')
parser.add_argument('--formulation_id', type=int, default=16, help='配方ID，默认为16')
parser.add_argument('--output_dir', type=str, help='输出目录，默认为test_output/formulation_<id>')

args = parser.parse_args()

# 设置Django环境
sys.path.append('/Users/xiaoji/Documents/AiSuan')  # 项目根目录
sys.path.append('/Users/xiaoji/Documents/AiSuan/backend')  # 后端目录
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'electrolyte_backend.settings')
django.setup()

# 导入需要的模块
from electrolyte.models import ElectrolyteFormulation, InputFile
from molyte_cursor.src.core.main import generate_input_files
from molyte_cursor.src.utils.logger import Logger

# 初始化日志
logger = Logger().get_logger()
logger.info("开始处理INP文件")

# 配方ID
formulation_id = args.formulation_id

try:
    # 获取配方和INP文件
    formulation = ElectrolyteFormulation.objects.get(id=formulation_id)
    input_file = InputFile.objects.get(formulation=formulation)
    
    print(f'找到配方: {formulation.name} (ID: {formulation_id}), INP文件路径: {input_file.file_path}')
    
    # 设置输出目录
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = f'/Users/xiaoji/Documents/AiSuan/test_output/formulation_{formulation_id}'
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f'准备使用molyte处理INP文件，输出目录: {output_dir}')
    
    # 读取INP文件内容
    with open(input_file.file_path, 'r', encoding='utf-8') as f:
        inp_content = f.read()
    
    # 准备配置
    config = {'inp_content': inp_content}
    
    # 使用molyte处理INP文件
    result = generate_input_files(config, output_dir)
    
    print(f'处理结果: {result}')
    print(f'输出目录: {output_dir}')
    
    # 列出生成的文件
    print('生成的文件:')
    for file_name in os.listdir(output_dir):
        file_path = os.path.join(output_dir, file_name)
        print(f'  - {file_name} ({os.path.getsize(file_path)} 字节)')
    
    # 显示文件内容摘要
    for file_name in ['input.inp', 'input.lammps']:
        file_path = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            print(f'\n{file_name} 文件内容预览:')
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                # 显示前10行
                for i, line in enumerate(lines[:10]):
                    print(f'  {i+1}: {line.rstrip()}')
                if len(lines) > 10:
                    print(f'  ... (总共 {len(lines)} 行)')

except Exception as e:
    print(f'处理过程中出错: {str(e)}')
    import traceback
    traceback.print_exc() 