"""
可视化模块

提供各种分析结果的可视化功能，包括RDF、MSD、溶剂分布等图表生成
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from typing import Dict, List, Tuple, Optional, Union, Any
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

class Visualizer:
    """
    可视化器
    
    提供各种分析结果的可视化功能，生成不同类型的图表
    """
    def __init__(self, output_dir: str, dpi: int = 300, logger=None):
        """
        初始化可视化器

        Args:
            output_dir: 输出图像的目录
            dpi: 图像分辨率
            logger: 日志记录器实例
        """
        self.output_dir = output_dir
        self.dpi = dpi
        self.logger = logger or logging.getLogger(__name__)
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 设置中文字体
        self._setup_fonts()
    
    def _setup_fonts(self):
        """
        设置Matplotlib字体，支持中文显示
        """
        try:
            # 尝试设置中文字体
            plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans', 'sans-serif']
            plt.rcParams['axes.unicode_minus'] = False  # 正确显示负号
            self.logger.debug("成功设置中文字体")
        except Exception as e:
            self.logger.warning(f"设置中文字体时出错: {e}")
    
    def plot_rdf(self, rdf_file: str, rdf_labels: List[str], output_filename: Optional[str] = None) -> str:
        """
        绘制径向分布函数(RDF)图像

        Args:
            rdf_file: RDF数据文件路径
            rdf_labels: RDF和CN标签列表
            output_filename: 输出文件名，如果为None则自动生成

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制RDF图像: {rdf_file}")
        
        # 设置输出文件名
        if output_filename is None:
            basename = os.path.basename(rdf_file).split('.')[0]
            output_filename = f"rdf_{basename}.png"
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            with open(rdf_file, 'r') as f:
                line3 = ''
                title3 = []
                current_rdf = []

                # 读取前4行注释
                for i in range(0, 4):
                    if i == 2:
                        line3 = f.readline()  # 第3行包含标题
                    else:
                        f.readline()  # 跳过其他行

                self.logger.debug(f"读取的标题行: {line3}")
                for idx, value in enumerate(line3.split()):
                    current_rdf.append([])
                    title3.append(value)

                # 读取数据、提取数据存入数组
                for line in f.readlines():
                    for idx, value in enumerate(line.split()):
                        current_rdf[idx].append(float(value))

            # 确保rdf_labels和current_rdf的数量匹配
            if len(rdf_labels) < len(current_rdf) // 2:
                self.logger.warning(f"rdf_labels太短: {len(rdf_labels)}，current_rdf数据数量: {len(current_rdf)}")
                return ""

            # 绘制RDF图
            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(111)
            self.logger.debug(f"标签数量: {len(rdf_labels)}")

            # 使用包含分子名的RDF标签进行绘制
            for idx in range(1, len(rdf_labels) // 2 + 1):
                rdf_label = rdf_labels[2 * (idx - 1)]  # 取出RDF标签
                self.logger.debug(f"正在绘制: {rdf_label}")
                ax.plot(current_rdf[1], current_rdf[2 * idx], linewidth=1.5, marker='_', 
                      markersize=3, color=f'C{idx-1}', label=rdf_label)

            # 绘制CN图像
            ax2 = ax.twinx()
            for idx in range(1, len(rdf_labels) // 2 + 1):
                cn_label = rdf_labels[2 * (idx - 1) + 1]  # 取出CN标签
                self.logger.debug(f"正在绘制: {cn_label}")
                ax2.plot(current_rdf[1], current_rdf[2 * idx + 1], linewidth=1.5, linestyle='--',
                       marker='.', markersize=3, color=f'C{idx-1}', label=cn_label)

            # 设置轴和标签
            ax.set_xlabel('$r$($\AA$)')
            ax.set_xlim(0, 8)
            ax.set_ylim(0)
            ax.set_ylabel('RDF(r)')
            ax.legend(loc='upper left')

            ax2.set_ylabel('CN')
            ax2.set_ylim(0)
            ax2.legend(loc='upper right')

            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"RDF图像已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制RDF图像时出错: {e}")
            return ""
    
    def plot_msd(self, msd_file: str, output_filename: Optional[str] = None) -> str:
        """
        绘制均方位移(MSD)图像

        Args:
            msd_file: MSD数据文件路径
            output_filename: 输出文件名，如果为None则自动生成

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制MSD图像: {msd_file}")
        
        # 设置输出文件名
        if output_filename is None:
            basename = os.path.basename(msd_file).split('.')[0]
            output_filename = f"msd_{basename}.png"
        
        output_path = os.path.join(self.output_dir, output_filename)
        
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
            ax.set_ylabel('MSD ($\AA^2$)')
            ax.legend()
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # 添加标题
            plt.title('均方位移随时间的变化')

            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"MSD图像已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制MSD图像时出错: {e}")
            return ""
    
    def plot_combined_msd(self, msd_files: Dict[str, str], output_filename: str, 
                          legend_labels: Optional[Dict[str, str]] = None) -> str:
        """
        将多个MSD文件的数据组合到一个图表中

        Args:
            msd_files: 样品名称到MSD文件路径的映射
            output_filename: 输出文件名
            legend_labels: 图例标签映射，如果为None则使用样品名称

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制组合MSD图像: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
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
            ax.set_ylabel('MSD ($\AA^2$)')
            ax.legend(loc='best')
            ax.grid(True, linestyle='--', alpha=0.7)
            
            # 添加标题
            plt.title('不同样品的均方位移(MSD)对比')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"组合MSD图像已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制组合MSD图像时出错: {e}")
            return ""
    
    def plot_diffusion_coefficients(self, diffusion_data: Dict[str, Dict[str, float]], 
                                   output_filename: str, sort_by: str = 'D_total') -> str:
        """
        绘制多个样品的扩散系数对比图

        Args:
            diffusion_data: 样品名称到扩散系数的映射，每个样品包含D_x, D_y, D_z, D_total
            output_filename: 输出文件名
            sort_by: 用于排序的扩散系数类型

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制扩散系数对比图: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 将数据转换为DataFrame
            data = []
            for sample, coeffs in diffusion_data.items():
                for direction, value in coeffs.items():
                    data.append({
                        '样品': sample,
                        '方向': direction.replace('D_', ''),
                        '扩散系数': value
                    })
            
            df = pd.DataFrame(data)
            
            # 排序样品
            sample_order = []
            if sort_by:
                # 根据指定的扩散系数排序
                sort_direction = sort_by.replace('D_', '')
                temp_df = df[df['方向'] == sort_direction].sort_values('扩散系数', ascending=False)
                sample_order = temp_df['样品'].tolist()
            
            # 绘制柱状图
            plt.figure(figsize=(12, 8))
            
            # 使用seaborn进行绘图
            ax = sns.barplot(data=df, x='样品', y='扩散系数', hue='方向', 
                            order=sample_order if sample_order else None)
            
            # 设置标签和标题
            plt.xlabel('样品')
            plt.ylabel('扩散系数 (10^-9 m^2/s)')
            plt.title('不同样品的扩散系数对比')
            
            # 旋转x轴标签，避免重叠
            plt.xticks(rotation=45, ha='right')
            
            # 添加数值标签
            for p in ax.patches:
                ax.annotate(f"{p.get_height():.2f}", 
                           (p.get_x() + p.get_width() / 2., p.get_height()),
                           ha='center', va='center', 
                           xytext=(0, 5), 
                           textcoords='offset points')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"扩散系数对比图已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制扩散系数对比图时出错: {e}")
            return ""
    
    def plot_solvent_distribution(self, solvent_data: Dict[str, Any], 
                                 output_filename: str, plot_type: str = 'bar') -> str:
        """
        绘制溶剂分布图

        Args:
            solvent_data: 溶剂分布数据
            output_filename: 输出文件名
            plot_type: 图表类型，可选'bar'、'pie'或'heatmap'

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制溶剂分布图: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            if plot_type == 'pie':
                # 饼图展示总体溶剂分布
                plt.figure(figsize=(10, 10))
                
                # 合并所有阳离子的平均溶剂化数据
                all_solvents = {}
                for cation_id, solvents in solvent_data['avg_cation_solvation'].items():
                    for solvent, count in solvents.items():
                        all_solvents[solvent] = all_solvents.get(solvent, 0) + count
                
                # 绘制饼图
                labels = all_solvents.keys()
                sizes = all_solvents.values()
                plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
                plt.axis('equal')
                plt.title('阳离子周围溶剂分布')
                
            elif plot_type == 'heatmap':
                # 热图展示离子和溶剂的关系
                plt.figure(figsize=(12, 8))
                
                # 构建热图数据
                cations = list(solvent_data['avg_cation_solvation'].keys())
                
                # 找出所有出现的溶剂类型
                all_solvent_types = set()
                for solvents in solvent_data['avg_cation_solvation'].values():
                    all_solvent_types.update(solvents.keys())
                all_solvent_types = sorted(all_solvent_types)
                
                # 创建热图数据
                heatmap_data = []
                for cation in cations:
                    row = []
                    for solvent in all_solvent_types:
                        row.append(solvent_data['avg_cation_solvation'][cation].get(solvent, 0))
                    heatmap_data.append(row)
                
                # 绘制热图
                sns.heatmap(heatmap_data, annot=True, fmt=".2f", 
                          xticklabels=all_solvent_types, yticklabels=cations)
                plt.title('阳离子-溶剂分布热图')
                plt.xlabel('溶剂类型')
                plt.ylabel('阳离子ID')
                
            else:  # 默认为柱状图
                plt.figure(figsize=(12, 8))
                
                # 构建柱状图数据
                data = []
                for cation_id, solvents in solvent_data['avg_cation_solvation'].items():
                    for solvent, count in solvents.items():
                        data.append({
                            '阳离子': f"Cation_{cation_id}",
                            '溶剂类型': solvent,
                            '平均溶剂化数': count
                        })
                
                df = pd.DataFrame(data)
                
                # 绘制柱状图
                ax = sns.barplot(data=df, x='阳离子', y='平均溶剂化数', hue='溶剂类型')
                
                # 添加数值标签
                for p in ax.patches:
                    ax.annotate(f"{p.get_height():.2f}", 
                               (p.get_x() + p.get_width() / 2., p.get_height()),
                               ha='center', va='center', 
                               xytext=(0, 5), 
                               textcoords='offset points')
                
                plt.title('阳离子周围溶剂分布')
                plt.xlabel('阳离子')
                plt.ylabel('平均溶剂化数')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"溶剂分布图已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制溶剂分布图时出错: {e}")
            return ""
    
    def plot_gaussian_energies(self, energy_data: Dict[str, float], output_filename: str) -> str:
        """
        绘制高斯计算能量对比图

        Args:
            energy_data: 分子名称到能量的映射
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制高斯能量对比图: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 排序数据
            sorted_data = sorted(energy_data.items(), key=lambda x: x[1])
            molecules, energies = zip(*sorted_data)
            
            # 能量转换为相对能量（相对于最低值）
            min_energy = min(energies)
            relative_energies = [e - min_energy for e in energies]
            
            # 绘制柱状图
            plt.figure(figsize=(12, 8))
            bars = plt.bar(molecules, relative_energies, color=plt.cm.tab10.colors)
            
            # 添加数值标签
            for i, bar in enumerate(bars):
                plt.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01,
                    f"{relative_energies[i]:.2f}",
                    ha='center', va='bottom'
                )
            
            # 设置标签和标题
            plt.xlabel('分子')
            plt.ylabel('相对能量 (Hartree)')
            plt.title('分子能量对比')
            
            # 旋转x轴标签，避免重叠
            plt.xticks(rotation=45, ha='right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"高斯能量对比图已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制高斯能量对比图时出错: {e}")
            return ""
    
    def plot_thermochemistry(self, thermo_data: Dict[str, Dict[str, float]], output_filename: str) -> str:
        """
        绘制热力学数据对比图

        Args:
            thermo_data: 分子名称到热力学数据的映射
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制热力学数据对比图: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 选择关键的热力学数据
            key_properties = [
                'zero_point_energy', 
                'thermal_enthalpy', 
                'thermal_free_energy'
            ]
            
            # 将数据转换为DataFrame
            data = []
            for molecule, properties in thermo_data.items():
                for prop in key_properties:
                    if prop in properties:
                        data.append({
                            '分子': molecule,
                            '热力学属性': prop,
                            '值': properties[prop]
                        })
            
            df = pd.DataFrame(data)
            
            # 绘制热力学数据对比
            plt.figure(figsize=(12, 8))
            
            # 使用seaborn进行绘图
            sns.barplot(data=df, x='分子', y='值', hue='热力学属性')
            
            # 设置标签和标题
            plt.xlabel('分子')
            plt.ylabel('能量 (Hartree)')
            plt.title('分子热力学属性对比')
            
            # 旋转x轴标签，避免重叠
            plt.xticks(rotation=45, ha='right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"热力学数据对比图已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制热力学数据对比图时出错: {e}")
            return ""
    
    def create_comparison_dashboard(self, data: Dict[str, Any], output_filename: str) -> str:
        """
        创建综合对比面板，包含多种分析结果

        Args:
            data: 包含各种分析结果的字典
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"创建综合对比面板: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 设置大图布局
            fig = plt.figure(figsize=(20, 16))
            # 设置网格布局
            grid = plt.GridSpec(2, 2, wspace=0.3, hspace=0.3)
            
            # 如果有扩散系数数据
            if 'diffusion' in data:
                ax1 = fig.add_subplot(grid[0, 0])
                
                # 准备扩散系数数据
                diff_data = []
                for sample, coeffs in data['diffusion'].items():
                    for direction, value in coeffs.items():
                        diff_data.append({
                            '样品': sample,
                            '方向': direction.replace('D_', ''),
                            '扩散系数': value
                        })
                
                df_diff = pd.DataFrame(diff_data)
                
                # 绘制扩散系数
                sns.barplot(data=df_diff, x='样品', y='扩散系数', hue='方向', ax=ax1)
                ax1.set_title('扩散系数对比')
                ax1.set_xlabel('样品')
                ax1.set_ylabel('扩散系数 (10^-9 m^2/s)')
                ax1.tick_params(axis='x', rotation=45)
            
            # 如果有溶剂分布数据
            if 'solvation' in data:
                ax2 = fig.add_subplot(grid[0, 1])
                
                # 合并所有阳离子的平均溶剂化数据
                all_solvents = {}
                for cation_id, solvents in data['solvation']['avg_cation_solvation'].items():
                    for solvent, count in solvents.items():
                        all_solvents[solvent] = all_solvents.get(solvent, 0) + count
                
                # 绘制饼图
                labels = all_solvents.keys()
                sizes = all_solvents.values()
                ax2.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
                ax2.axis('equal')
                ax2.set_title('溶剂分布')
            
            # 如果有Gaussian能量数据
            if 'energies' in data:
                ax3 = fig.add_subplot(grid[1, 0])
                
                # 准备能量数据
                sorted_energy_data = sorted(data['energies'].items(), key=lambda x: x[1])
                molecules, energies = zip(*sorted_energy_data)
                
                # 能量转换为相对能量（相对于最低值）
                min_energy = min(energies)
                relative_energies = [(e - min_energy) * 627.5095 for e in energies]  # 转换为kcal/mol
                
                # 绘制柱状图
                bars = ax3.bar(molecules, relative_energies, color=plt.cm.tab10.colors)
                
                # 添加数值标签
                for i, bar in enumerate(bars):
                    ax3.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.1,
                        f"{relative_energies[i]:.2f}",
                        ha='center', va='bottom'
                    )
                
                ax3.set_title('相对能量 (kcal/mol)')
                ax3.set_xlabel('分子')
                ax3.set_ylabel('相对能量')
                ax3.tick_params(axis='x', rotation=45)
            
            # 如果有RDF数据
            if 'rdf_data' in data:
                ax4 = fig.add_subplot(grid[1, 1])
                
                # 假设是一个预先绘制好的RDF图的文件路径
                if os.path.exists(data['rdf_data']):
                    img = plt.imread(data['rdf_data'])
                    ax4.imshow(img)
                    ax4.axis('off')
                    ax4.set_title('RDF分析')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"综合对比面板已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"创建综合对比面板时出错: {e}")
            return ""
    
    def plot_ion_pair_analysis(self, ion_pair_data: Dict[str, Any], output_filename: str) -> str:
        """
        绘制离子对分析图表

        Args:
            ion_pair_data: 离子对分析数据，包含CIP、SSIP、AGG和自由阳离子的计数和百分比
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制离子对分析图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            
            # 提取数据
            pair_types = ['CIP', 'SSIP', 'AGG', '自由阳离子']
            counts = [
                ion_pair_data.get('cip_count', 0),
                ion_pair_data.get('ssip_count', 0),
                ion_pair_data.get('agg_count', 0),
                ion_pair_data.get('free_cation_count', 0)
            ]
            percentages = [
                ion_pair_data.get('cip_percentage', 0),
                ion_pair_data.get('ssip_percentage', 0),
                ion_pair_data.get('agg_percentage', 0),
                ion_pair_data.get('free_cation_percentage', 0)
            ]
            
            # 绘制计数柱状图
            ax1.bar(pair_types, counts, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
            ax1.set_title('离子对类型计数')
            ax1.set_xlabel('离子对类型')
            ax1.set_ylabel('计数')
            
            # 添加数值标签
            for i, count in enumerate(counts):
                ax1.text(i, count + 0.1, str(count), ha='center')
            
            # 绘制百分比饼图
            wedges, texts, autotexts = ax2.pie(
                percentages, 
                labels=pair_types,
                autopct='%1.1f%%',
                startangle=90,
                colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
            )
            ax2.set_title('离子对类型百分比')
            ax2.axis('equal')
            
            # 设置自动文本的字体大小
            plt.setp(autotexts, size=10, weight="bold")
            
            # 添加图例
            plt.figlegend(
                wedges,
                pair_types,
                loc="lower center",
                ncol=4,
                frameon=False,
            )
            
            # 添加总标题
            plt.suptitle('离子对分析结果', fontsize=16)
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"离子对分析图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制离子对分析图表时出错: {e}")
            return ""
    
    def plot_conductivity_data(self, electrochemical_data: Dict[str, Any], output_filename: str) -> str:
        """
        绘制电导率数据图表

        Args:
            electrochemical_data: 电化学数据
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制电导率数据图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, ax1 = plt.figure(figsize=(10, 7)), plt.gca()
            
            # 提取数据
            conductivity = electrochemical_data.get('conductivity_mS_cm', 0)
            molar_conductivity = electrochemical_data.get('molar_conductivity', 0)
            concentration = electrochemical_data.get('concentration', 1.0)
            sample_name = electrochemical_data.get('sample_name', 'Unknown')
            
            # 创建条形图
            bar_width = 0.35
            index = np.arange(1)
            
            # 电导率条形图
            ax1.bar(index, [conductivity], bar_width, label='电导率', color='#1f77b4')
            ax1.set_ylabel('电导率 (mS/cm)')
            ax1.set_title(f'样品 {sample_name} 电导率 (浓度 {concentration} mol/L)')
            ax1.set_xticks(index)
            ax1.set_xticklabels([sample_name])
            
            # 摩尔电导率条形图
            ax2 = ax1.twinx()
            ax2.bar(index + bar_width, [molar_conductivity], bar_width, label='摩尔电导率', color='#ff7f0e')
            ax2.set_ylabel('摩尔电导率 (S·m²/mol)')
            
            # 添加数值标签
            ax1.text(index[0], conductivity + 0.1, f"{conductivity:.3f}", ha='center')
            ax2.text(index[0] + bar_width, molar_conductivity + 0.1, f"{molar_conductivity:.3e}", ha='center')
            
            # 添加图例
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"电导率数据图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制电导率数据图表时出错: {e}")
            return ""
    
    def plot_viscosity_data(self, electrochemical_data: Dict[str, Any], output_filename: str) -> str:
        """
        绘制粘度数据图表

        Args:
            electrochemical_data: 电化学数据
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制粘度数据图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, ax1 = plt.figure(figsize=(9, 6)), plt.gca()
            
            # 提取数据
            viscosity = electrochemical_data.get('viscosity_mPa_s', 0)
            diffusion = electrochemical_data.get('diffusion_coeffs', {}).get('D_total', 0)
            sample_name = electrochemical_data.get('sample_name', 'Unknown')
            
            # 创建条形图
            bar_width = 0.35
            index = np.arange(1)
            
            # 粘度条形图
            ax1.bar(index, [viscosity], bar_width, label='粘度', color='#2ca02c')
            ax1.set_ylabel('粘度 (mPa·s)')
            ax1.set_title(f'样品 {sample_name} 粘度与扩散系数')
            ax1.set_xticks(index)
            ax1.set_xticklabels([sample_name])
            
            # 扩散系数条形图
            ax2 = ax1.twinx()
            ax2.bar(index + bar_width, [diffusion], bar_width, label='扩散系数', color='#d62728')
            ax2.set_ylabel('扩散系数 (10⁻⁹ m²/s)')
            
            # 添加数值标签
            ax1.text(index[0], viscosity + 0.1, f"{viscosity:.3f}", ha='center')
            ax2.text(index[0] + bar_width, diffusion + 0.1, f"{diffusion:.3f}", ha='center')
            
            # 添加图例
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"粘度数据图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制粘度数据图表时出错: {e}")
            return ""
    
    def plot_comparative_conductivity(self, electrochemical_data: Dict[str, Dict[str, Any]], 
                                    output_filename: str) -> str:
        """
        绘制比较性电导率数据图表

        Args:
            electrochemical_data: 多个样品的电化学数据
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制比较性电导率数据图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, ax1 = plt.figure(figsize=(12, 8)), plt.gca()
            
            # 提取数据
            samples = []
            conductivity_values = []
            molar_conductivity_values = []
            
            for sample, data in electrochemical_data.items():
                samples.append(sample)
                conductivity_values.append(data.get('conductivity_mS_cm', 0))
                molar_conductivity_values.append(data.get('molar_conductivity', 0))
            
            # 创建条形图
            bar_width = 0.35
            index = np.arange(len(samples))
            
            # 电导率条形图
            bars1 = ax1.bar(index, conductivity_values, bar_width, label='电导率', color='#1f77b4')
            ax1.set_ylabel('电导率 (mS/cm)')
            ax1.set_title('不同样品电导率比较')
            ax1.set_xticks(index + bar_width / 2)
            ax1.set_xticklabels(samples, rotation=45, ha='right')
            
            # 摩尔电导率条形图
            ax2 = ax1.twinx()
            bars2 = ax2.bar(index + bar_width, molar_conductivity_values, bar_width, label='摩尔电导率', color='#ff7f0e')
            ax2.set_ylabel('摩尔电导率 (S·m²/mol)')
            
            # 添加数值标签
            for i, (bar, value) in enumerate(zip(bars1, conductivity_values)):
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f"{value:.3f}", ha='center', va='bottom')
            
            for i, (bar, value) in enumerate(zip(bars2, molar_conductivity_values)):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f"{value:.3e}", ha='center', va='bottom')
            
            # 添加图例
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"比较性电导率数据图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制比较性电导率数据图表时出错: {e}")
            return ""
    
    def plot_comparative_viscosity(self, electrochemical_data: Dict[str, Dict[str, Any]], 
                                 output_filename: str) -> str:
        """
        绘制比较性粘度数据图表

        Args:
            electrochemical_data: 多个样品的电化学数据
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制比较性粘度数据图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, ax1 = plt.figure(figsize=(12, 8)), plt.gca()
            
            # 提取数据
            samples = []
            viscosity_values = []
            diffusion_values = []
            
            for sample, data in electrochemical_data.items():
                samples.append(sample)
                viscosity_values.append(data.get('viscosity_mPa_s', 0))
                diffusion_values.append(data.get('diffusion_coefficients', {}).get('D_total', 0))
            
            # 创建条形图
            bar_width = 0.35
            index = np.arange(len(samples))
            
            # 粘度条形图
            bars1 = ax1.bar(index, viscosity_values, bar_width, label='粘度', color='#2ca02c')
            ax1.set_ylabel('粘度 (mPa·s)')
            ax1.set_title('不同样品粘度比较')
            ax1.set_xticks(index + bar_width / 2)
            ax1.set_xticklabels(samples, rotation=45, ha='right')
            
            # 扩散系数条形图
            ax2 = ax1.twinx()
            bars2 = ax2.bar(index + bar_width, diffusion_values, bar_width, label='扩散系数', color='#d62728')
            ax2.set_ylabel('扩散系数 (10⁻⁹ m²/s)')
            
            # 添加数值标签
            for i, (bar, value) in enumerate(zip(bars1, viscosity_values)):
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f"{value:.3f}", ha='center', va='bottom')
            
            for i, (bar, value) in enumerate(zip(bars2, diffusion_values)):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f"{value:.3f}", ha='center', va='bottom')
            
            # 添加图例
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
            
            # 保存图像
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"比较性粘度数据图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制比较性粘度数据图表时出错: {e}")
            return ""
    
    def plot_comparative_transference(self, electrochemical_data: Dict[str, Dict[str, Any]], 
                                    output_filename: str) -> str:
        """
        绘制比较性离子迁移数图表

        Args:
            electrochemical_data: 多个样品的电化学数据
            output_filename: 输出文件名

        Returns:
            生成图像的路径
        """
        self.logger.info(f"绘制比较性离子迁移数图表: {output_filename}")
        
        output_path = os.path.join(self.output_dir, output_filename)
        
        try:
            # 创建图表
            fig, ax = plt.figure(figsize=(10, 7)), plt.gca()
            
            # 提取数据
            samples = []
            transference_values = []
            
            for sample, data in electrochemical_data.items():
                samples.append(sample)
                transference_values.append(data.get('transference_number', 0))
            
            # 创建条形图
            bars = ax.bar(samples, transference_values, color='#9467bd')
            ax.set_ylabel('阳离子迁移数')
            ax.set_title('不同样品阳离子迁移数比较')
            ax.set_ylim(0, 1)  # 迁移数范围在0-1之间
            ax.axhline(y=0.5, linestyle='--', color='gray', alpha=0.7)  # 添加0.5的参考线
            
            # 添加数值标签
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                       f"{height:.3f}", ha='center', va='bottom')
            
            # 保存图像
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(output_path, dpi=self.dpi)
            self.logger.info(f"比较性离子迁移数图表已保存到: {output_path}")
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"绘制比较性离子迁移数图表时出错: {e}")
            return "" 