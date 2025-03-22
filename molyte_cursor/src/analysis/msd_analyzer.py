"""
MSD分析模块

用于分析Mean Square Displacement (MSD)和相关计算
"""

import os
import logging
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple, Optional, Union, Any
import re

class MSDAnalyzer:
    """
    均方位移(MSD)分析器
    
    提供处理和分析分子动力学模拟中均方位移的功能
    """
    def __init__(self, work_dir: str, logger=None):
        """
        初始化MSD分析器

        Args:
            work_dir: 工作目录路径
            logger: 日志记录器实例
        """
        self.working_dir = work_dir
        self.logger = logger or logging.getLogger(__name__)
    
    def plot_msd(self, msd_file: str, output_filename: str) -> None:
        """
        读取MSD文件，并绘制MSD曲线图

        Args:
            msd_file: MSD数据文件路径
            output_filename: 输出图像文件路径
        """
        self.logger.info(f"绘制MSD图像: {msd_file} -> {output_filename}")
        
        try:
            with open(msd_file, 'r') as f:
                line1 = ''
                title1 = []
                current_msd = []

                # 读取前3行注释
                for i in range(0, 3):
                    if i == 0:
                        line1 = f.readline()  # 第1行包含标题
                    else:
                        f.readline()  # 跳过其他行

                self.logger.debug(f"读取的标题行: {line1}")
                for idx, value in enumerate(line1.split()):
                    current_msd.append([])
                    if idx == 1:
                        title1.append(f'{value}_X')
                    elif idx == 2:
                        title1.append(f'{value}_Y')
                    elif idx == 3:
                        title1.append(f'{value}_Z')
                    else:
                        title1.append(value)

                # 读取数据、提取数据存入数组
                for line in f.readlines():
                    for idx, value in enumerate(line.split()):
                        current_msd[idx].append(float(value))

            # 将时间单位转换为ps
            time_in_ps = [t / 1000 for t in current_msd[0]]  # 假设时间单位是fs

            # 绘制MSD图
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111)

            # 分别绘制X、Y、Z方向和总体MSD
            ax.plot(time_in_ps, current_msd[1], 'r-', label='X方向')
            ax.plot(time_in_ps, current_msd[2], 'g-', label='Y方向')
            ax.plot(time_in_ps, current_msd[3], 'b-', label='Z方向')
            if len(current_msd) > 4:  # 如果有总体MSD数据
                ax.plot(time_in_ps, current_msd[4], 'k-', linewidth=2, label='总体')

            # 设置轴和标签
            ax.set_xlabel('时间 (ps)')
            ax.set_ylabel(r'MSD ($\AA^2$)')
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.7)

            # 保存图像
            plt.tight_layout()
            plt.savefig(output_filename, dpi=300)
            self.logger.info(f"MSD图像已保存到: {output_filename}")
            plt.close()
            
        except Exception as e:
            self.logger.error(f"绘制MSD图像时出错: {e}")
            raise
    
    def calculate_diffusion_coefficient(self, msd_file: str, time_unit: str = 'fs', skip_initial: int = 10) -> Dict[str, float]:
        """
        计算扩散系数
        
        使用Einstein关系从MSD数据计算扩散系数：D = MSD/(6*t)
        
        Args:
            msd_file: MSD数据文件路径
            time_unit: 时间单位，默认为'fs'
            skip_initial: 跳过初始的数据点数量，因为初始阶段可能不是线性的
            
        Returns:
            各方向和总体的扩散系数（单位：10^-9 m^2/s）
        """
        self.logger.info(f"计算扩散系数: {msd_file}")
        
        try:
            # 读取MSD数据
            with open(msd_file, 'r') as f:
                line1 = ''
                current_msd = []

                # 读取前3行注释
                for i in range(0, 3):
                    if i == 0:
                        line1 = f.readline()  # 第1行包含标题
                    else:
                        f.readline()  # 跳过其他行

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
                time_factor = 1e-15  # fs到s的转换
            elif time_unit == 'ps':
                time_factor = 1e-12  # ps到s的转换
            elif time_unit == 'ns':
                time_factor = 1e-9   # ns到s的转换
            
            # 跳过初始数据点
            times = np.array(current_msd[0][skip_initial:]) * time_factor
            msd_x = np.array(current_msd[1][skip_initial:])
            msd_y = np.array(current_msd[2][skip_initial:])
            msd_z = np.array(current_msd[3][skip_initial:])
            msd_total = np.array(current_msd[4][skip_initial:]) if len(current_msd) > 4 else msd_x + msd_y + msd_z
            
            # 线性拟合
            slope_x, _ = np.polyfit(times, msd_x, 1)
            slope_y, _ = np.polyfit(times, msd_y, 1)
            slope_z, _ = np.polyfit(times, msd_z, 1)
            slope_total, _ = np.polyfit(times, msd_total, 1)
            
            # 计算扩散系数 D = slope/6 (适用于3D)
            # 转换为常用单位：10^-9 m^2/s
            d_x = slope_x / 6 * 1e9
            d_y = slope_y / 6 * 1e9
            d_z = slope_z / 6 * 1e9
            d_total = slope_total / 6 * 1e9
            
            results = {
                'D_x': d_x,
                'D_y': d_y,
                'D_z': d_z,
                'D_total': d_total
            }
            
            self.logger.info(f"扩散系数计算结果: {results}")
            return results
            
        except Exception as e:
            self.logger.error(f"计算扩散系数时出错: {e}")
            raise
    
    def plot_combined_msd(self, msd_files: Dict[str, str], output_filename: str, legend_labels: Optional[Dict[str, str]] = None) -> None:
        """
        将多个MSD文件的数据组合到一个图表中
        
        Args:
            msd_files: 样品名称到MSD文件路径的映射
            output_filename: 输出图像文件路径
            legend_labels: 图例标签映射，如果为None则使用样品名称
        """
        self.logger.info(f"绘制组合MSD图像 -> {output_filename}")
        
        try:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111)
            
            colors = plt.cm.tab10.colors  # 使用Matplotlib的颜色循环
            markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
            
            for idx, (sample_name, msd_file) in enumerate(msd_files.items()):
                if not os.path.exists(msd_file):
                    self.logger.warning(f"MSD文件不存在: {msd_file}")
                    continue
                
                # 读取MSD数据
                with open(msd_file, 'r') as f:
                    current_msd = [[] for _ in range(5)]  # 假设最多有5列数据
                    
                    # 跳过前3行注释
                    for _ in range(3):
                        f.readline()
                    
                    # 读取数据
                    for line in f.readlines():
                        parts = line.split()
                        for i in range(min(len(parts), 5)):
                            current_msd[i].append(float(parts[i]))
                
                # 将时间单位转换为ps
                time_in_ps = [t / 1000 for t in current_msd[0]]  # 假设时间单位是fs
                
                # 使用总体MSD（通常是第5列）
                if len(current_msd) >= 5 and len(current_msd[4]) > 0:
                    label = legend_labels.get(sample_name, sample_name) if legend_labels else sample_name
                    color = colors[idx % len(colors)]
                    marker = markers[idx % len(markers)]
                    
                    # 绘制数据点（仅绘制部分点，避免过度拥挤）
                    step = max(1, len(time_in_ps) // 20)
                    ax.plot(time_in_ps[::step], current_msd[4][::step], 
                           marker=marker, markersize=5, linestyle='-', 
                           color=color, label=label, markevery=step)
            
            # 设置轴和标签
            ax.set_xlabel('时间 (ps)')
            ax.set_ylabel(r'MSD ($\AA^2$)')
            ax.legend(loc='best')
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # 添加标题
            plt.title('不同样品的均方位移(MSD)对比')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_filename, dpi=300)
            self.logger.info(f"组合MSD图像已保存到: {output_filename}")
            plt.close()
            
        except Exception as e:
            self.logger.error(f"绘制组合MSD图像时出错: {e}")
            raise
            
    def process_msd_data(self, sample_name: str, current_directory: str, output_directory: str) -> Dict[str, float]:
        """
        处理单个样品的MSD数据，生成图像并计算扩散系数
        
        Args:
            sample_name: 样品名称
            current_directory: 当前处理目录
            output_directory: 输出数据目录
            
        Returns:
            计算得到的扩散系数
        """
        self.logger.info(f"处理样品 {sample_name} 的MSD数据")
        
        try:
            # 定义文件路径
            msd_file = os.path.join(current_directory, 'out_msd.dat')
            msd_image_path = os.path.join(output_directory, f'msd_{sample_name}.png')
            diffusion_file = os.path.join(output_directory, f'diffusion_{sample_name}.txt')
            
            # 检查MSD文件是否存在
            if not os.path.exists(msd_file):
                self.logger.warning(f"MSD文件不存在: {msd_file}")
                return {}
            
            # 绘制MSD图像
            self.plot_msd(msd_file, msd_image_path)
            
            # 计算扩散系数
            diffusion_coeffs = self.calculate_diffusion_coefficient(msd_file)
            
            # 保存扩散系数到文件
            with open(diffusion_file, 'w') as f:
                f.write(f"样品: {sample_name}\n")
                f.write("扩散系数 (10^-9 m^2/s):\n")
                f.write(f"X方向: {diffusion_coeffs['D_x']:.6f}\n")
                f.write(f"Y方向: {diffusion_coeffs['D_y']:.6f}\n")
                f.write(f"Z方向: {diffusion_coeffs['D_z']:.6f}\n")
                f.write(f"总体: {diffusion_coeffs['D_total']:.6f}\n")
            
            self.logger.info(f"样品 {sample_name} 的MSD数据处理完成")
            return diffusion_coeffs
            
        except Exception as e:
            self.logger.error(f"处理样品 {sample_name} 的MSD数据时出错: {e}")
            raise
            
    def calculate_ionic_conductivity(self, msd_file: str, temperature: float = 298.15, 
                                   concentration: float = 1.0, 
                                   cation_valence: int = 1, anion_valence: int = 1) -> Dict[str, float]:
        """
        使用Nernst-Einstein方程计算离子电导率
        
        电导率 σ = (e^2 * N_A * C / kB / T) * (z_+^2 * D_+ + z_-^2 * D_-)
        
        Args:
            msd_file: MSD数据文件路径
            temperature: 温度 (K)
            concentration: 电解质浓度 (mol/L)
            cation_valence: 阳离子价数
            anion_valence: 阴离子价数
            
        Returns:
            包含电导率和相关参数的字典
        """
        self.logger.info(f"计算离子电导率: {msd_file}")
        
        try:
            # 计算扩散系数
            diffusion_coeffs = self.calculate_diffusion_coefficient(msd_file)
            
            # 物理常数
            e_charge = 1.60217662e-19  # 电子电荷 (C)
            kb = 1.38064852e-23  # 玻尔兹曼常数 (J/K)
            na = 6.022e23  # 阿伏伽德罗常数 (mol^-1)
            
            # 转换浓度从mol/L到mol/m^3
            concentration_mol_m3 = concentration * 1000  # 1000 L/m^3
            
            # 转换扩散系数从10^-9 m^2/s到m^2/s
            d_cation = diffusion_coeffs['D_total'] * 1e-9  # 使用总体扩散系数作为阳离子扩散系数
            d_anion = diffusion_coeffs['D_total'] * 1e-9   # 使用总体扩散系数作为阴离子扩散系数
            
            # 计算电导率 (S/m)
            conductivity = (e_charge**2 * na * concentration_mol_m3 / (kb * temperature)) * \
                          (cation_valence**2 * d_cation + anion_valence**2 * d_anion)
            
            # 计算摩尔电导率 (S*m^2/mol)
            molar_conductivity = conductivity / concentration_mol_m3
            
            # 计算阳离子迁移数
            cation_transport = cation_valence**2 * d_cation
            anion_transport = anion_valence**2 * d_anion
            transference_number = cation_transport / (cation_transport + anion_transport)
            
            results = {
                'conductivity': conductivity,  # S/m
                'conductivity_mS_cm': conductivity * 0.1,  # mS/cm
                'molar_conductivity': molar_conductivity,  # S*m^2/mol
                'transference_number': transference_number,  # 无量纲
                'diffusion_cation': d_cation,  # m^2/s
                'diffusion_anion': d_anion,  # m^2/s
                'temperature': temperature,  # K
                'concentration': concentration  # mol/L
            }
            
            self.logger.info(f"离子电导率计算结果: {results['conductivity_mS_cm']:.4f} mS/cm")
            return results
            
        except Exception as e:
            self.logger.error(f"计算离子电导率时出错: {e}")
            raise
            
    def process_electrochemical_data(self, sample_name: str, current_directory: str, output_directory: str,
                                   temperature: float = 298.15, concentration: float = 1.0,
                                   cation_valence: int = 1, anion_valence: int = 1,
                                   particle_radius: float = 0.2) -> Dict[str, Any]:
        """
        处理样品的电化学数据，计算各种电化学性能参数
        
        Args:
            sample_name: 样品名称
            current_directory: 当前处理目录
            output_directory: 输出数据目录
            temperature: 温度 (K)
            concentration: 电解质浓度 (mol/L)
            cation_valence: 阳离子价数
            anion_valence: 阴离子价数
            particle_radius: 粒子半径 (nm)，用于估算粘度
            
        Returns:
            包含所有电化学性能参数的字典
        """
        self.logger.info(f"处理样品 {sample_name} 的电化学数据")
        
        try:
            # 定义文件路径
            msd_file = os.path.join(current_directory, 'out_msd.dat')
            electrochemical_file = os.path.join(output_directory, f'electrochemical_{sample_name}.txt')
            
            # 检查MSD文件是否存在
            if not os.path.exists(msd_file):
                self.logger.warning(f"MSD文件不存在: {msd_file}")
                return {}
            
            # 计算扩散系数
            diffusion_coeffs = self.calculate_diffusion_coefficient(msd_file)
            
            # 计算电导率和相关参数
            conductivity_results = self.calculate_ionic_conductivity(
                msd_file, temperature, concentration, cation_valence, anion_valence
            )
            
            # 物理常数
            kb = 1.38064852e-23  # 玻尔兹曼常数 (J/K)
            pi = 3.14159265359  # 圆周率
            
            # 转换半径从nm到m
            radius_m = particle_radius * 1e-9
            
            # 转换扩散系数从10^-9 m^2/s到m^2/s
            d_total = diffusion_coeffs['D_total'] * 1e-9
            
            # 使用Stokes-Einstein关系计算粘度 (Pa*s)
            viscosity = kb * temperature / (6 * pi * radius_m * d_total)
            
            # 整合所有结果
            results = {
                'sample_name': sample_name,
                'diffusion_coeffs': diffusion_coeffs,  # 10^-9 m^2/s
                'conductivity': conductivity_results['conductivity'],  # S/m
                'conductivity_mS_cm': conductivity_results['conductivity_mS_cm'],  # mS/cm
                'molar_conductivity': conductivity_results['molar_conductivity'],  # S*m^2/mol
                'transference_number': conductivity_results['transference_number'],  # 无量纲
                'viscosity': viscosity,  # Pa*s
                'viscosity_mPa_s': viscosity * 1000,  # mPa*s (cP)
                'temperature': temperature,  # K
                'concentration': concentration  # mol/L
            }
            
            # 保存结果到文件
            with open(electrochemical_file, 'w') as f:
                f.write(f"样品: {sample_name}\n")
                f.write(f"温度: {temperature} K\n")
                f.write(f"浓度: {concentration} mol/L\n\n")
                
                f.write("扩散系数 (10^-9 m^2/s):\n")
                f.write(f"X方向: {diffusion_coeffs['D_x']:.6f}\n")
                f.write(f"Y方向: {diffusion_coeffs['D_y']:.6f}\n")
                f.write(f"Z方向: {diffusion_coeffs['D_z']:.6f}\n")
                f.write(f"总体: {diffusion_coeffs['D_total']:.6f}\n\n")
                
                f.write("电化学性能参数:\n")
                f.write(f"电导率: {results['conductivity']:.6e} S/m ({results['conductivity_mS_cm']:.4f} mS/cm)\n")
                f.write(f"摩尔电导率: {results['molar_conductivity']:.6e} S*m^2/mol\n")
                f.write(f"阳离子迁移数: {results['transference_number']:.4f}\n")
                f.write(f"粘度: {results['viscosity']:.6e} Pa*s ({results['viscosity_mPa_s']:.4f} mPa*s)\n")
            
            self.logger.info(f"样品 {sample_name} 的电化学数据处理完成")
            return results
            
        except Exception as e:
            self.logger.error(f"处理样品 {sample_name} 的电化学数据时出错: {e}")
            raise 