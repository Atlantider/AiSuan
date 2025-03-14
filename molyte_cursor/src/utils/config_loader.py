"""
配置加载器模块，负责加载和管理配置信息
"""
import os
import yaml
from pathlib import Path

class ConfigLoader:
    """配置加载器类"""
    
    def __init__(self, config_dir=None):
        """初始化配置加载器
        
        Args:
            config_dir: 配置目录路径，如果为None则使用默认路径
        """
        if config_dir is None:
            # 使用当前文件所在包的上级目录中的config文件夹
            self.config_dir = Path(__file__).parents[3] / "config"
        else:
            self.config_dir = Path(config_dir)
        
        self.paths_config = None
        self.simulation_params = None
        
    def load_paths_config(self, config_file="default_paths.yaml"):
        """加载路径配置
        
        Args:
            config_file: 配置文件名
            
        Returns:
            加载的路径配置字典
        """
        config_path = self.config_dir / "paths" / config_file
        with open(config_path, "r") as f:
            self.paths_config = yaml.safe_load(f)
        return self.paths_config
    
    def load_simulation_params(self, config_file="default_params.yaml"):
        """加载模拟参数配置
        
        Args:
            config_file: 配置文件名
            
        Returns:
            加载的模拟参数配置字典
        """
        config_path = self.config_dir / "simulation_params" / config_file
        with open(config_path, "r") as f:
            self.simulation_params = yaml.safe_load(f)
        return self.simulation_params
    
    def get_path(self, path_key):
        """获取特定路径
        
        Args:
            path_key: 路径键名
            
        Returns:
            对应的路径对象
        """
        if self.paths_config is None:
            self.load_paths_config()
        
        path_str = self.paths_config["paths"].get(path_key)
        if path_str:
            return Path(path_str)
        return None
    
    def get_file_path(self, file_key):
        """获取特定文件路径
        
        Args:
            file_key: 文件键名
            
        Returns:
            对应的文件路径对象
        """
        if self.paths_config is None:
            self.load_paths_config()
        
        file_str = self.paths_config["files"].get(file_key)
        if file_str:
            return Path(file_str)
        return None
    
    def get_simulation_param(self, param_path):
        """获取特定模拟参数
        
        Args:
            param_path: 参数路径，以点分隔，如 'npt.steps'
            
        Returns:
            对应的参数值
        """
        if self.simulation_params is None:
            self.load_simulation_params()
        
        # 按点分割参数路径
        path_parts = param_path.split('.')
        
        # 从simulation字典开始导航
        current = self.simulation_params.get("simulation", {})
        
        # 逐级导航
        for part in path_parts:
            if part in current:
                current = current[part]
            else:
                return None
                
        return current 