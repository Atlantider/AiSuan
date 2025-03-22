"""
Slurm作业监控工具

提供用于监控和管理Slurm作业的实用函数
"""
import os
import re
import time
import subprocess
import logging
from typing import Dict, List, Optional, Tuple, Any

logger = logging.getLogger(__name__)

def check_job_status(job_id: str) -> str:
    """
    检查Slurm作业的状态
    
    Args:
        job_id: Slurm作业ID
    
    Returns:
        作业状态 (PENDING, RUNNING, COMPLETED, FAILED, CANCELLED, 或 UNKNOWN)
    """
    try:
        cmd = f"sacct -j {job_id} -o State -n"
        result = subprocess.run(cmd, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if result.returncode != 0:
            logger.error(f"无法获取作业状态: {result.stderr.decode().strip()}")
            return "UNKNOWN"
        
        output = result.stdout.decode().strip()
        
        # 提取状态
        if "PENDING" in output:
            return "PENDING"
        elif "RUNNING" in output:
            return "RUNNING"
        elif "COMPLETED" in output:
            return "COMPLETED"
        elif "FAILED" in output or "TIMEOUT" in output:
            return "FAILED"
        elif "CANCELLED" in output:
            return "CANCELLED"
        else:
            return "UNKNOWN"
    except Exception as e:
        logger.error(f"检查作业状态时出错: {str(e)}")
        return "UNKNOWN"

def wait_for_job_completion(job_id: str, polling_interval: int = 60, max_wait_time: Optional[int] = None) -> Tuple[bool, str]:
    """
    等待Slurm作业完成
    
    Args:
        job_id: Slurm作业ID
        polling_interval: 轮询间隔（秒）
        max_wait_time: 最大等待时间（秒），如果为None则无限等待
        
    Returns:
        (是否成功, 最终状态)
    """
    start_time = time.time()
    logger.info(f"开始等待作业 {job_id} 完成")
    
    while True:
        status = check_job_status(job_id)
        
        if status in ["COMPLETED", "FAILED", "CANCELLED"]:
            logger.info(f"作业 {job_id} 已完成，状态: {status}")
            return status == "COMPLETED", status
        
        # 检查是否超过最大等待时间
        if max_wait_time and (time.time() - start_time > max_wait_time):
            logger.warning(f"等待作业 {job_id} 超时")
            return False, "TIMEOUT"
        
        logger.debug(f"作业 {job_id} 状态: {status}，继续等待")
        time.sleep(polling_interval)

def cancel_job(job_id: str) -> bool:
    """
    取消Slurm作业
    
    Args:
        job_id: Slurm作业ID
        
    Returns:
        是否成功取消
    """
    try:
        cmd = f"scancel {job_id}"
        result = subprocess.run(cmd, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if result.returncode != 0:
            logger.error(f"取消作业失败: {result.stderr.decode().strip()}")
            return False
        
        logger.info(f"已取消作业 {job_id}")
        return True
    except Exception as e:
        logger.error(f"取消作业时出错: {str(e)}")
        return False

def get_job_info(job_id: str) -> Dict[str, Any]:
    """
    获取Slurm作业的详细信息
    
    Args:
        job_id: Slurm作业ID
        
    Returns:
        包含作业信息的字典
    """
    try:
        cmd = f"sacct -j {job_id} -o JobID,JobName,State,Elapsed,CPUTime,NodeList,ReqMem -P"
        result = subprocess.run(cmd, shell=True, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        if result.returncode != 0:
            logger.error(f"获取作业信息失败: {result.stderr.decode().strip()}")
            return {"error": result.stderr.decode().strip()}
        
        output = result.stdout.decode().strip()
        lines = output.split('\n')
        
        if len(lines) < 2:
            return {"error": "无法解析作业信息"}
        
        # 解析表头和数据
        headers = lines[0].split('|')
        values = lines[1].split('|')
        
        job_info = {}
        for i, header in enumerate(headers):
            if i < len(values):
                job_info[header.strip()] = values[i].strip()
        
        return job_info
    except Exception as e:
        logger.error(f"获取作业信息时出错: {str(e)}")
        return {"error": str(e)}

def get_job_log(job_id: str, log_type: str = "out") -> str:
    """
    获取Slurm作业的日志内容
    
    Args:
        job_id: Slurm作业ID
        log_type: 日志类型 ('out' 或 'err')
        
    Returns:
        日志内容
    """
    try:
        # 首先获取作业信息以确定日志文件位置
        job_info = get_job_info(job_id)
        
        # 典型的Slurm日志文件命名格式
        log_pattern = f"slurm-{job_id}.{log_type}"
        
        # 尝试在几个常见位置查找日志文件
        possible_locations = [
            os.getcwd(),
            os.path.expanduser("~"),
            "/tmp",
        ]
        
        for location in possible_locations:
            log_path = os.path.join(location, log_pattern)
            if os.path.exists(log_path):
                with open(log_path, 'r', encoding='utf-8', errors='replace') as f:
                    return f.read()
        
        return f"找不到作业 {job_id} 的{log_type}日志文件"
    except Exception as e:
        logger.error(f"获取作业日志时出错: {str(e)}")
        return f"获取日志时出错: {str(e)}" 