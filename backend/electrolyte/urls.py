from django.urls import path, include
from rest_framework.routers import DefaultRouter
from rest_framework.decorators import api_view
from rest_framework.response import Response
from .views import (
    SolventViewSet,
    SaltViewSet,
    ElectrolyteFormulationViewSet,
    CalculationViewSet,
    save_input_file,
    check_input_file,
    get_calculation_status,
    get_calculation_results,
    submit_calculation,
    process_web_recipe,
    fix_formulation_smile,
    fix_all_smiles
)
from .test_views import test_file_generation, test_smile_fix

# 添加测试API视图
@api_view(['GET'])
def test_api(request):
    return Response({'status': 'ok', 'message': 'API服务正常'})

# 创建路由器并注册视图集
router = DefaultRouter()
router.register(r'solvents', SolventViewSet)
router.register(r'salts', SaltViewSet)
router.register(r'formulations', ElectrolyteFormulationViewSet, basename='formulation')
router.register(r'calculations', CalculationViewSet, basename='calculation')

urlpatterns = [
    path('', include(router.urls)),
    
    # 测试API路由
    path('test-api/', test_api, name='test-api'),
    
    # 新增INP文件相关路由 - 修改为接受任何ID类型
    path('formulations/<str:formulation_id>/input-file/', save_input_file, name='save-input-file'),
    path('formulations/<str:formulation_id>/check-input-file/', check_input_file, name='check-input-file'),
    
    # 新增计算任务相关路由
    path('submit-calculation/', submit_calculation, name='submit-calculation'),
    path('calculations/<str:calculation_id>/status/', get_calculation_status, name='calculation-status'),
    path('calculations/<str:calculation_id>/results/', get_calculation_results, name='calculation-results'),
    
    # 新增从网页直接提交配方的路由
    path('process-recipe/', process_web_recipe, name='process-web-recipe'),
    
    # 修复SMILE字符串的路由
    path('formulations/<str:pk>/fix-smile/', fix_formulation_smile, name='fix-formulation-smile'),
    path('fix-all-smiles/', fix_all_smiles, name='fix-all-smiles'),
    
    # 测试路由
    path('test/file-generation/', test_file_generation, name='test_file_generation'),
    path('test/smile-fix/', test_smile_fix, name='test_smile_fix'),
]
