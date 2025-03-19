"""
RDF分析模块

用于分析Radius Distribution Function (RDF)和相关计算
"""

import os
import logging
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple, Optional, Union

class RDFAnalyzer:
    """
    径向分布函数(RDF)分析器
    
    提供处理和分析分子动力学模拟中径向分布函数的功能
    """
    def __init__(self, working_dir: str, logger=None):
        """
        初始化RDF分析器

        Args:
            working_dir: 工作目录路径
            logger: 日志记录器实例
        """
        self.working_dir = working_dir
        self.logger = logger or logging.getLogger(__name__)
    
    def assign_atoms_to_molecules(self, atom_counts: Dict[str, int]) -> Dict[str, Tuple[int, int]]:
        """
        根据每个分子的原子数量，计算分子中原子的索引范围

        Args:
            atom_counts: 分子名称到原子数量的映射

        Returns:
            分子名称到(start_idx, end_idx)的映射
        """
        atom_ranges = {}
        start_idx = 1  # 从1开始，因为索引从1开始

        for molecule, count in atom_counts.items():
            end_idx = start_idx + count - 1
            atom_ranges[molecule] = (start_idx, end_idx)
            self.logger.debug(f"指定分子 '{molecule}' 到原子索引 {start_idx}-{end_idx}")
            start_idx = end_idx + 1

        self.logger.debug(f"最终原子范围: {atom_ranges}")
        return atom_ranges

    def find_molecule_for_atom(self, atom_index: int, atom_ranges: Dict[str, Tuple[int, int]]) -> Optional[str]:
        """
        根据原子索引，查找该原子属于哪个分子

        Args:
            atom_index: 原子的类型编号
            atom_ranges: 分子名称到(start_idx, end_idx)的映射

        Returns:
            分子名称或None（如果未找到）
        """
        for molecule, (start, end) in atom_ranges.items():
            if start <= atom_index <= end:
                self.logger.debug(f"原子类型 {atom_index} 属于分子 '{molecule}'")
                return molecule
        self.logger.debug(f"原子类型 {atom_index} 不属于任何分子")
        return None

    def parse_in_list_with_molecule(self, in_list_path: str, atom_counts: Dict[str, int]) -> Tuple[List[str], List[str]]:
        """
        解析指定的输入列表文件，提取元素列表和RDF对，并为原子添加所属分子标签

        Args:
            in_list_path: 输入列表文件的路径
            atom_counts: 分子名称到原子数量的映射

        Returns:
            包含元素列表和RDF标签的元组
        """
        self.logger.info(f"解析输入列表文件: {in_list_path}")
        
        # 分配原子到分子
        atom_ranges = self.assign_atoms_to_molecules(atom_counts)
        
        element_lists = []
        rdf_labels = []
        
        try:
            with open(in_list_path, 'r') as f:
                lines = f.readlines()
                
                for line in lines:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split()
                    
                    # 处理元素列表
                    if len(parts) >= 2 and parts[0] == "element":
                        element_lists.append(line)
                        self.logger.debug(f"找到元素列表: {line}")
                    
                    # 处理RDF
                    elif len(parts) >= 2 and parts[0] == "rdf":
                        # 解析原子类型号
                        type1 = int(parts[1])
                        type2 = int(parts[2])
                        
                        # 获取原子所属的分子
                        mol1 = self.find_molecule_for_atom(type1, atom_ranges)
                        mol2 = self.find_molecule_for_atom(type2, atom_ranges)
                        
                        # 构建带有分子信息的标签
                        label = f"RDF_{parts[1]}_{parts[2]}"
                        cn_label = f"CN_{parts[1]}_{parts[2]}"
                        
                        if mol1 and mol2:
                            label = f"RDF_{mol1}({parts[1]})_{mol2}({parts[2]})"
                            cn_label = f"CN_{mol1}({parts[1]})_{mol2}({parts[2]})"
                        
                        self.logger.debug(f"生成RDF标签: {label}, CN标签: {cn_label}")
                        rdf_labels.append(label)
                        rdf_labels.append(cn_label)
            
            self.logger.info(f"解析完成，找到 {len(element_lists)} 个元素列表和 {len(rdf_labels)//2} 个RDF对")
            return element_lists, rdf_labels
            
        except Exception as e:
            self.logger.error(f"解析输入列表文件时出错: {e}")
            raise

    def replace_rdf_labels_in_file(self, original_file: str, rdf_labels: List[str], output_file: str) -> None:
        """
        替换文件中的RDF标签并生成新的输出文件

        Args:
            original_file: 原始RDF文件路径
            rdf_labels: RDF和CN标签列表
            output_file: 输出文件路径
        """
        self.logger.info(f"替换文件中的RDF标签: {original_file} -> {output_file}")
        
        try:
            # 读取原始文件
            with open(original_file, 'r') as f:
                lines = f.readlines()
            
            # 确保有足够的行
            if len(lines) < 4:
                self.logger.error(f"原始文件格式不正确，行数不足: {len(lines)}")
                return
            
            # 替换第三行的标签（索引为2）
            header_parts = lines[2].strip().split()
            modified_header = []
            
            # 保留第一和第二列（通常是步数和距离）
            if len(header_parts) >= 2:
                modified_header.extend(header_parts[:2])
            
            # 添加自定义RDF标签
            modified_header.extend(rdf_labels)
            
            # 更新标题行
            lines[2] = " ".join(modified_header) + "\n"
            
            # 写入输出文件
            with open(output_file, 'w') as f:
                f.writelines(lines)
                
            self.logger.info(f"标签替换完成，已保存到: {output_file}")
            
        except Exception as e:
            self.logger.error(f"替换RDF标签时出错: {e}")
            raise

    def plot_rdf(self, rdf_file: str, rdf_labels: List[str], output_filename: str) -> None:
        """
        读取RDF文件，并根据传入的rdf_labels绘制RDF和CN图像

        Args:
            rdf_file: RDF数据文件路径
            rdf_labels: RDF和CN标签列表
            output_filename: 输出图像文件路径
        """
        self.logger.info(f"绘制RDF图像: {rdf_file} -> {output_filename}")
        
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
                return

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
            ax.set_xlabel(r'$r$($\AA$)')
            ax.set_xlim(0, 8)
            ax.set_ylim(0)
            ax.set_ylabel('RDF(r)')
            ax.legend(loc='upper left')

            ax2.set_ylabel('CN')
            ax2.set_ylim(0)
            ax2.legend(loc='upper right')

            # 保存图像
            plt.tight_layout()
            plt.savefig(output_filename, dpi=300)
            self.logger.info(f"RDF图像已保存到: {output_filename}")
            plt.close()
            
        except Exception as e:
            self.logger.error(f"绘制RDF图像时出错: {e}")
            raise 