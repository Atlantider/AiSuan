from django.urls import path, include
from rest_framework.routers import DefaultRouter
from .views import (
    SolventViewSet,
    SaltViewSet,
    ElectrolyteFormulationViewSet,
    CalculationViewSet,
    save_input_file,
    check_input_file,
    get_calculation_status,
    get_calculation_results,
    submit_calculation
)
from .test_views import test_file_generation

# 创建路由器并注册视图集
router = DefaultRouter()
router.register(r'solvents', SolventViewSet)
router.register(r'salts', SaltViewSet)
router.register(r'formulations', ElectrolyteFormulationViewSet, basename='formulation')
router.register(r'calculations', CalculationViewSet, basename='calculation')

urlpatterns = [
    path('', include(router.urls)),
    
    # 新增INP文件相关路由 - 修改为接受任何ID类型
    path('formulations/<str:formulation_id>/input-file/', save_input_file, name='save-input-file'),
    path('formulations/<str:formulation_id>/check-input-file/', check_input_file, name='check-input-file'),
    
    # 新增计算任务相关路由
    path('submit-calculation/', submit_calculation, name='submit-calculation'),
    path('calculations/<str:calculation_id>/status/', get_calculation_status, name='calculation-status'),
    path('calculations/<str:calculation_id>/results/', get_calculation_results, name='calculation-results'),
    
    # 测试路由
    path('test/file-generation/', test_file_generation, name='test_file_generation'),
]
