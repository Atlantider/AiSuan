"""
命令执行器模块，用于执行系统命令并处理结果
"""
import subprocess
from pathlib import Path
from .logger import Logger

class CommandExecutor:
    """命令执行器类"""
    
    def __init__(self):
        """初始化命令执行器"""
        self.logger = Logger().get_logger()
    
    def run_command(self, command, cwd=None, check=True, capture_output=True, shell=True):
        """运行shell命令，并记录输出和错误信息
        
        Args:
            command: 要执行的命令
            cwd: 工作目录
            check: 是否检查命令执行状态
            capture_output: 是否捕获输出
            shell: 是否使用shell执行命令
            
        Returns:
            执行结果对象
            
        Raises:
            subprocess.CalledProcessError: 当命令执行失败且check=True时
        """
        self.logger.info(f"Running command: {command}")
        
        try:
            result = subprocess.run(
                command, 
                shell=shell, 
                check=check, 
                stdout=subprocess.PIPE if capture_output else None, 
                stderr=subprocess.PIPE if capture_output else None,
                cwd=cwd
            )
            
            if capture_output:
                output = result.stdout.decode().strip()
                if output:
                    self.logger.info(output)
                    
            return result
            
        except subprocess.CalledProcessError as e:
            if capture_output:
                error_msg = e.stderr.decode().strip()
                self.logger.error(f"Command failed: {error_msg}")
            raise
    
    def run_command_safe(self, command, cwd=None, capture_output=True, shell=True):
        """安全地运行shell命令，不会抛出异常
        
        Args:
            command: 要执行的命令
            cwd: 工作目录
            capture_output: 是否捕获输出
            shell: 是否使用shell执行命令
            
        Returns:
            (成功标志, 执行结果对象)
        """
        try:
            result = self.run_command(command, cwd, True, capture_output, shell)
            return True, result
        except subprocess.CalledProcessError as e:
            return False, e 