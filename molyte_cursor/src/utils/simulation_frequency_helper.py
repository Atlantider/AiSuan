"""
自适应模拟输出频率计算工具

该模块提供用于计算分子动力学模拟中自适应输出频率的函数。
根据模拟步长、平衡步数、生产步数和系统大小，动态调整输出频率。
"""

import math
import logging

# 设置日志
logger = logging.getLogger("SimulationFrequencyHelper")

def calculate_adaptive_frequencies(formulation_data, total_particles=None):
    """
    计算自适应输出频率
    
    根据模拟参数和系统大小动态计算热力学数据和轨迹输出的频率。
    
    Args:
        formulation_data (dict): 包含模拟参数的配方数据字典
        total_particles (int, optional): 系统中的总粒子数。如果为None，将尝试从formulation_data估计
        
    Returns:
        tuple: (thermo_freq, trj_freq_npt, trj_freq_nvt) 热力学数据频率和两个阶段的轨迹输出频率
    """
    parameters = formulation_data.get('parameters', {})
    
    # 获取模拟参数(使用默认值作为后备)
    time_step = float(parameters.get('time_step', 1.0))  # fs
    equilibration_steps = int(parameters.get('equilibration_steps', 10000))
    production_steps = int(parameters.get('production_steps', 100000))
    
    # 系统大小(粒子数)估计
    if total_particles is None:
        total_particles = 0
        for component_type in ['cations', 'anions', 'solvents']:
            for component in formulation_data.get(component_type, []):
                # 尝试获取数量，支持'number'和'count'两种键名
                count = component.get('number', component.get('count', 0))
                try:
                    total_particles += int(count)
                except (TypeError, ValueError):
                    # 如果转换失败，记录警告并继续
                    logger.warning(f"无法将组分数量转换为整数: {count}")
    
    logger.info(f"估计的系统总粒子数: {total_particles}")
    
    # 目标数据点数量 - 根据时间步长略微调整
    time_step_factor = math.sqrt(time_step)  # 时间步长的平方根作为调整因子
    
    target_thermo_points = 200 * time_step_factor  # 热力学数据期望点数
    target_equil_frames = 50 * time_step_factor    # 平衡阶段期望帧数
    target_prod_frames = 100 * time_step_factor    # 生产阶段期望帧数
    
    # 系统大小调整因子(大系统应该输出更少帧)
    size_factor = max(1.0, min(3.0, total_particles / 1000))
    
    # 计算自适应频率(向上取整到10的倍数)
    thermo_freq = math.ceil((equilibration_steps + production_steps) / target_thermo_points / 10) * 10
    trj_freq_npt = math.ceil(equilibration_steps / target_equil_frames / 10) * 10 * size_factor
    trj_freq_nvt = math.ceil(production_steps / target_prod_frames / 10) * 10 * size_factor
    
    # 确保频率在合理范围内
    thermo_freq = max(100, min(5000, int(thermo_freq)))
    trj_freq_npt = max(500, min(10000, int(trj_freq_npt)))
    trj_freq_nvt = max(1000, min(20000, int(trj_freq_nvt)))
    
    logger.info(f"自适应输出频率: thermo_freq={thermo_freq}, trj_freq_npt={trj_freq_npt}, trj_freq_nvt={trj_freq_nvt}")
    logger.info(f"基于: time_step={time_step}fs, equil_steps={equilibration_steps}, prod_steps={production_steps}")
    
    return thermo_freq, trj_freq_npt, trj_freq_nvt

def get_standard_frequencies():
    """
    获取标准固定频率(用作后备)
    
    Returns:
        tuple: (thermo_freq, trj_freq_npt, trj_freq_nvt) 标准输出频率
    """
    return 1000, 2000, 10000 