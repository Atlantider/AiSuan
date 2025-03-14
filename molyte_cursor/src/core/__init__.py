"""
核心模块，包含主程序入口和集成工作流
"""

from .main import main
from .integrated_workflow import IntegratedWorkflow

__all__ = [
    'main',
    'IntegratedWorkflow'
]
