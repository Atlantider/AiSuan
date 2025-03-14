"""
日志工具模块，提供统一的日志记录功能
"""
import logging
import sys
from pathlib import Path

class Logger:
    """日志管理类"""
    
    _instance = None
    
    def __new__(cls, *args, **kwargs):
        """单例模式实现"""
        if cls._instance is None:
            cls._instance = super(Logger, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self, log_level=logging.INFO, log_file=None):
        """初始化日志管理器
        
        Args:
            log_level: 日志级别
            log_file: 日志文件路径，如果为None则只输出到控制台
        """
        if self._initialized:
            return
            
        self.logger = logging.getLogger('molyte')
        self.logger.setLevel(log_level)
        self.logger.propagate = False
        
        # 清除已有的处理器
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
        
        # 创建控制台处理器
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(console_handler)
        
        # 如果指定了日志文件，创建文件处理器
        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(levelname)s - %(message)s'
            ))
            self.logger.addHandler(file_handler)
        
        self._initialized = True
    
    def get_logger(self):
        """获取日志记录器
        
        Returns:
            日志记录器对象
        """
        return self.logger
        
    def info(self, message):
        """记录信息级别日志
        
        Args:
            message: 日志消息
        """
        self.logger.info(message)
        
    def warning(self, message):
        """记录警告级别日志
        
        Args:
            message: 日志消息
        """
        self.logger.warning(message)
        
    def error(self, message):
        """记录错误级别日志
        
        Args:
            message: 日志消息
        """
        self.logger.error(message)
        
    def debug(self, message):
        """记录调试级别日志
        
        Args:
            message: 日志消息
        """
        self.logger.debug(message) 