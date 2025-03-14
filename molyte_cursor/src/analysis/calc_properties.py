#!/usr/bin/env python
"""
计算和可视化分子电子性质的命令行工具
"""
import argparse
import sys
import os
from pathlib import Path
import numpy as np
from scipy import stats
import pandas as pd

# 添加父目录到sys.path，使模块导入正常工作
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

from molyte_cursor.src.analysis.gaussian_processor import GaussianProcessor
from molyte_cursor.src.utils.logger import Logger

# 物理常数
KB = 1.38064852e-23  # 玻尔兹曼常数 (J/K)
E_CHARGE = 1.60217662e-19  # 电子电荷 (C)
EPSILON_0 = 8.85418782e-12  # 真空介电常数 (F/m)
NA = 6.022e23  # 阿伏伽德罗常数 (mol^-1)
PI = np.pi  # 圆周率

def calculate_ionic_conductivity(msd_file: str, temperature: float = 298.15, 
                                time_unit: str = 'fs', concentration: float = 1.0,
                                cation_valence: int = 1, anion_valence: int = 1,
                                skip_initial: int = 10) -> float:
    """
    使用Nernst-Einstein方程计算离子电导率
    
    电导率 σ = (e^2 * N_A * C / kB / T) * (z_+^2 * D_+ + z_-^2 * D_-)
    
    Args:
        msd_file: MSD数据文件路径
        temperature: 温度 (K)
        time_unit: 时间单位，可选 'fs', 'ps', 'ns'
        concentration: 电解质浓度 (mol/L)
        cation_valence: 阳离子价数
        anion_valence: 阴离子价数
        skip_initial: 跳过的初始数据点数
        
    Returns:
        离子电导率 (S/m)
    """
    # 读取MSD数据
    with open(msd_file, 'r') as f:
        line1 = ''
        current_msd = []

        # 读取前3行注释
        for i in range(0, 3):
            if i == 0:
                line1 = f.readline()
            else:
                f.readline()

        # 初始化数据数组
        for idx, _ in enumerate(line1.split()):
            current_msd.append([])

        # 读取数据
        for line in f.readlines():
            for idx, value in enumerate(line.split()):
                current_msd[idx].append(float(value))
    
    # 转换时间单位为秒
    time_factor = 1.0
    if time_unit == 'fs':
        time_factor = 1e-15
    elif time_unit == 'ps':
        time_factor = 1e-12
    elif time_unit == 'ns':
        time_factor = 1e-9
    
    # 跳过初始数据点
    times = np.array(current_msd[0][skip_initial:]) * time_factor
    msd_cation = np.array(current_msd[4][skip_initial:]) if len(current_msd) > 4 else None
    msd_anion = np.array(current_msd[4][skip_initial:]) if len(current_msd) > 4 else None
    
    # 线性拟合计算扩散系数
    if msd_cation is not None:
        slope_cation, _ = np.polyfit(times, msd_cation, 1)
        d_cation = slope_cation / 6  # 扩散系数 (m^2/s)
    else:
        # 如果没有分别的阳离子数据，使用总体MSD
        msd_total = np.array(current_msd[4][skip_initial:])
        slope_total, _ = np.polyfit(times, msd_total, 1)
        d_cation = slope_total / 6  # 估计值
    
    if msd_anion is not None:
        slope_anion, _ = np.polyfit(times, msd_anion, 1)
        d_anion = slope_anion / 6  # 扩散系数 (m^2/s)
    else:
        # 如果没有分别的阴离子数据，使用总体MSD或与阳离子相同
        d_anion = d_cation
    
    # 转换浓度从mol/L到mol/m^3
    concentration_mol_m3 = concentration * 1000  # 1000 L/m^3
    
    # 计算电导率 (S/m)
    conductivity = (E_CHARGE**2 * NA * concentration_mol_m3 / (KB * temperature)) * \
                  (cation_valence**2 * d_cation + anion_valence**2 * d_anion)
    
    return conductivity

def calculate_dielectric_constant(dipole_file: str, temperature: float = 298.15, 
                                volume: float = 1.0, dipole_unit: str = 'e*A') -> float:
    """
    使用系统偶极矩涨落计算介电常数
    
    介电常数 ε = 1 + <M^2>/(3*ε_0*V*k_B*T)
    
    Args:
        dipole_file: 偶极矩数据文件路径
        temperature: 温度 (K)
        volume: 模拟体系体积 (nm^3)
        dipole_unit: 偶极矩单位，可选 'e*A', 'Debye'
        
    Returns:
        介电常数 (无量纲)
    """
    # 读取偶极矩数据
    dipole_data = np.loadtxt(dipole_file, skiprows=1)
    
    # 提取偶极矩分量
    if dipole_data.ndim > 1 and dipole_data.shape[1] >= 4:
        # 假设文件格式：时间 Mx My Mz
        mx = dipole_data[:, 1]
        my = dipole_data[:, 2]
        mz = dipole_data[:, 3]
    else:
        raise ValueError("偶极矩数据格式错误，需要至少4列：时间 Mx My Mz")
    
    # 计算总偶极矩
    m_total = np.sqrt(mx**2 + my**2 + mz**2)
    
    # 单位转换
    dipole_factor = 1.0
    if dipole_unit == 'e*A':
        # 从e*A到C*m的转换
        dipole_factor = E_CHARGE * 1e-10
    elif dipole_unit == 'Debye':
        # 从Debye到C*m的转换（1 D = 3.33564e-30 C*m）
        dipole_factor = 3.33564e-30
    
    # 转换体积从nm^3到m^3
    volume_m3 = volume * 1e-27
    
    # 计算偶极矩平方的平均值
    m_squared_avg = np.mean(m_total**2) * (dipole_factor**2)
    
    # 计算介电常数
    dielectric_constant = 1 + m_squared_avg / (3 * EPSILON_0 * volume_m3 * KB * temperature)
    
    return dielectric_constant

def calculate_transference_number(cation_diffusion: float, anion_diffusion: float,
                                cation_valence: int = 1, anion_valence: int = 1) -> float:
    """
    计算离子迁移数
    
    迁移数 t_+ = (z_+*D_+)/(z_+*D_+ + z_-*D_-)
    
    Args:
        cation_diffusion: 阳离子扩散系数 (m^2/s)
        anion_diffusion: 阴离子扩散系数 (m^2/s)
        cation_valence: 阳离子价数
        anion_valence: 阴离子价数
        
    Returns:
        阳离子迁移数 (无量纲，范围0-1)
    """
    numerator = cation_valence * cation_diffusion
    denominator = cation_valence * cation_diffusion + anion_valence * anion_diffusion
    
    if denominator == 0:
        return 0.5  # 默认值
    
    transference_number = numerator / denominator
    return transference_number

def calculate_haven_ratio(conductivity_md: float, conductivity_ne: float) -> float:
    """
    计算Haven比率
    
    Haven比率 = σ_MD / σ_NE
    
    Args:
        conductivity_md: 从MD轨迹直接计算的电导率
        conductivity_ne: 使用Nernst-Einstein方程计算的电导率
        
    Returns:
        Haven比率
    """
    if conductivity_ne == 0:
        return 1.0  # 默认值
    
    haven_ratio = conductivity_md / conductivity_ne
    return haven_ratio

def calculate_diffusion_activation_energy(diffusion_coeffs: list, temperatures: list) -> tuple:
    """
    使用Arrhenius方程计算扩散活化能
    
    D = D0 * exp(-Ea/kB*T)
    
    Args:
        diffusion_coeffs: 不同温度下的扩散系数列表 (m^2/s)
        temperatures: 对应的温度列表 (K)
        
    Returns:
        活化能 (kJ/mol) 和指前因子 D0
    """
    # 转换为Arrhenius形式
    x = 1 / (np.array(temperatures) * KB / E_CHARGE)  # 1/kT in eV^-1
    y = np.log(np.array(diffusion_coeffs))
    
    # 线性拟合
    slope, intercept, r_value, _, _ = stats.linregress(x, y)
    
    # 计算活化能 (kJ/mol)
    activation_energy = -slope * E_CHARGE * NA / 1000
    
    # 计算指前因子
    d0 = np.exp(intercept)
    
    return activation_energy, d0, r_value**2

def calculate_viscosity(diffusion_coeff: float, particle_radius: float = 0.1, 
                       temperature: float = 298.15) -> float:
    """
    使用Stokes-Einstein关系计算粘度
    
    η = kB*T / (6π*r*D)
    
    Args:
        diffusion_coeff: 粒子扩散系数 (m^2/s)
        particle_radius: 粒子半径 (nm)
        temperature: 温度 (K)
        
    Returns:
        粘度 (Pa*s)
    """
    # 转换半径从nm到m
    radius_m = particle_radius * 1e-9
    
    # 使用Stokes-Einstein方程计算粘度
    viscosity = KB * temperature / (6 * PI * radius_m * diffusion_coeff)
    
    return viscosity

def main():
    """命令行入口函数"""
    # 配置参数解析器
    parser = argparse.ArgumentParser(
        description="计算和可视化分子电子性质（HOMO, LUMO, ESP等）"
    )
    
    # 添加参数
    parser.add_argument(
        "input",
        help="输入Excel文件路径，包含溶剂分子信息"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        default="gaussian_calc",
        help="输出目录，默认为'gaussian_calc'"
    )
    
    parser.add_argument(
        "-r", "--run",
        action="store_true",
        help="生成输入文件后立即运行高斯计算"
    )
    
    parser.add_argument(
        "-g", "--gaussian",
        default="g16",
        help="高斯可执行文件路径或命令，默认为'g16'"
    )
    
    parser.add_argument(
        "-m", "--method",
        default="B3LYP",
        help="计算方法，默认为'B3LYP'"
    )
    
    parser.add_argument(
        "-b", "--basis",
        default="6-31G(d)",
        help="基组，默认为'6-31G(d)'"
    )
    
    parser.add_argument(
        "-c", "--calc-type",
        default="opt freq",
        help="计算类型，默认为'opt freq'"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="输出详细日志信息"
    )
    
    parser.add_argument(
        "-p", "--parse-only",
        action="store_true",
        help="只解析已有的高斯输出文件而不生成新的输入"
    )
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置日志级别
    logger = Logger(log_level="DEBUG" if args.verbose else "INFO").get_logger()
    
    # 创建高斯处理器
    processor = GaussianProcessor()
    
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    
    # 检查输入路径是否存在
    if not input_path.exists():
        logger.error(f"输入路径不存在: {input_path}")
        sys.exit(1)
    
    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 根据模式执行不同操作
    if args.parse_only:
        # 只解析已有的输出文件
        logger.info(f"解析目录 {output_dir} 中的高斯输出文件")
        results, summary_df = processor.parse_all_outputs(output_dir)
        
        if results:
            # 可视化结果
            plot_files = processor.visualize_results(results, output_dir)
            
            if summary_df is not None:
                logger.info(f"解析了 {len(results)} 个文件，结果摘要已保存到 {output_dir}/gaussian_results_summary.csv")
                
                # 打印前5行摘要数据
                print("\n前5行摘要数据:")
                print(summary_df.head().to_string())
            
            if plot_files:
                logger.info(f"生成了 {len(plot_files)} 个可视化图表")
                for plot_file in plot_files:
                    logger.info(f"  - {plot_file}")
        else:
            logger.error(f"在目录 {output_dir} 中未找到任何高斯输出文件")
            sys.exit(1)
    else:
        # 完整流程：从Excel生成输入、可选运行计算、解析结果、可视化
        logger.info(f"开始从Excel文件 {input_path} 计算分子性质")
        
        # 执行计算和可视化流程
        summary_df, plot_files = processor.calculate_and_visualize_properties(
            input_path,
            output_dir,
            run_gaussian=args.run,
            gaussian_exe=args.gaussian,
            calc_type=args.calc_type,
            basis_set=args.basis,
            method=args.method
        )
        
        if summary_df is not None:
            logger.info("处理完成，结果摘要:")
            print("\n摘要数据:")
            print(summary_df.to_string())
        
        if plot_files:
            logger.info(f"生成了 {len(plot_files)} 个可视化图表:")
            for plot_file in plot_files:
                logger.info(f"  - {plot_file}")
        
        if args.run:
            logger.info(f"所有计算和分析已完成，结果保存在目录 {output_dir}")
        else:
            logger.info(f"已生成高斯输入文件和批处理脚本，请按照提示运行脚本进行计算")

if __name__ == "__main__":
    main() 