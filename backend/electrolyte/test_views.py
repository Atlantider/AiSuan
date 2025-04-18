"""
测试视图模块
"""
import os
import sys
import json
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
from rest_framework.decorators import api_view, permission_classes
from rest_framework.permissions import IsAuthenticated
import logging
from pathlib import Path
from .models import ElectrolyteFormulation
from .utils import fix_smile_in_formulation, fix_all_ec_smile_errors

logger = logging.getLogger(__name__)

# 打印当前系统路径
print("Django中的Python路径:")
for p in sys.path:
    print(f"  - {p}")

# 添加molyte_cursor的路径
project_root = '/Users/xiaoji/Documents/AiSuan'
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    print(f"已添加 {project_root} 到Python路径")

# 尝试导入 molyte_cursor 模块
try:
    from molyte_cursor.src.core.main import generate_input_files
    from molyte_cursor.src.io.file_generator import LAMMPSFileGenerator
    MOLYTE_AVAILABLE = True
    print("成功导入molyte_cursor模块!")
except ImportError as e:
    MOLYTE_AVAILABLE = False
    print(f"导入错误: {e}")

@csrf_exempt
def test_file_generation(request):
    """
    测试文件生成功能的视图函数
    
    POST参数:
        - config: 包含配方信息的JSON配置
    """
    if request.method != 'POST':
        return JsonResponse({'error': '只支持POST请求'}, status=405)
    
    # 检查molyte_cursor是否可用
    if not MOLYTE_AVAILABLE:
        return JsonResponse({'error': 'molyte_cursor模块不可用'}, status=500)
    
    try:
        # 解析请求体中的JSON数据
        try:
            data = json.loads(request.body)
            config = data.get('config')
            if not config:
                # 使用默认配置
                config = {
                    "formulation_name": "测试视图配方",
                    "temperature": 300,
                    "box_size": 40,
                    "concentration": 1.0,
                    "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
                    "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
                }
        except json.JSONDecodeError:
            # 使用默认配置
            config = {
                "formulation_name": "测试视图配方",
                "temperature": 300,
                "box_size": 40,
                "concentration": 1.0,
                "salts": [{"name": "LiPF6", "cation": "Li", "anion": "PF6", "concentration": 1.0}],
                "solvents": [{"name": "EC", "smile": "C1OC(=O)O1", "concentration": 1.0}]
            }
        
        # 创建输出目录
        output_dir = os.path.join(settings.MEDIA_ROOT, 'api_test_output')
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成文件
        result = generate_input_files(config, output_dir)
        
        # 计算相对于MEDIA_URL的文件路径
        media_url = settings.MEDIA_URL
        inp_relative_path = os.path.relpath(result['inp_file'], settings.MEDIA_ROOT)
        lammps_relative_path = os.path.relpath(result['lammps_file'], settings.MEDIA_ROOT)
        
        # 构建返回数据
        response_data = {
            'success': True,
            'message': '文件生成成功',
            'files': {
                'inp_file': f"{media_url}{inp_relative_path}",
                'lammps_file': f"{media_url}{lammps_relative_path}"
            },
            'config': config
        }
        
        return JsonResponse(response_data)
        
    except Exception as e:
        return JsonResponse({'error': f'文件生成失败: {str(e)}'}, status=500)

@api_view(['GET'])
@permission_classes([IsAuthenticated])
def test_smile_fix(request):
    """
    测试修复SMILE字符串功能
    """
    try:
        # 获取所有配方的数量
        total_count = ElectrolyteFormulation.objects.count()
        
        # 尝试修复所有配方
        fixed_count, error_count, details = fix_all_ec_smile_errors()
        
        # 构建测试报告
        result = {
            "success": True,
            "message": f"完成SMILE字符串测试修复",
            "total_formulations": total_count,
            "fixed_count": fixed_count,
            "error_count": error_count,
            "details": details[:10]  # 只返回前10个结果，避免响应过大
        }
        
        if len(details) > 10:
            result["note"] = f"结果被截断，共 {len(details)} 个配方已处理"
        
        return JsonResponse(result)
    except Exception as e:
        logger.error(f"测试SMILE修复功能时出错: {str(e)}")
        return JsonResponse({
            "success": False,
            "message": f"测试失败: {str(e)}"
        }, status=500) 