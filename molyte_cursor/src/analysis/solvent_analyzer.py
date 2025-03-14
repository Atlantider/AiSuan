"""
溶剂分析器模块，用于分析溶剂化结构、处理周期性边界条件并导出溶剂化结构
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
import re
import matplotlib.pyplot as plt
from ..utils.logger import Logger
from .trajectory_processor import TrajectoryProcessor

class SolventAnalyzer:
    """溶剂分析器类，用于分析溶剂化结构并处理周期性边界条件"""
    
    def __init__(self):
        """初始化溶剂分析器"""
        self.logger = Logger().get_logger()
        self.trajectory_processor = TrajectoryProcessor()
    
    def read_xyz_trajectory(self, trajectory_file, frame_index=-1):
        """读取XYZ轨迹文件中的指定帧
        
        Args:
            trajectory_file: 轨迹文件路径
            frame_index: 帧索引，-1表示最后一帧
            
        Returns:
            (元素列表, 坐标数组, 盒子尺寸)
        """
        self.logger.info(f"从文件 {trajectory_file} 读取第 {frame_index} 帧")
        
        frames = []
        elements = []
        current_frame = []
        atoms_per_frame = None
        box_size = None
        
        with open(trajectory_file, 'r') as f:
            lines = f.readlines()
            
            line_index = 0
            while line_index < len(lines):
                if not lines[line_index].strip():
                    line_index += 1
                    continue
                
                # 读取原子数
                atom_count = int(lines[line_index].strip())
                if atoms_per_frame is None:
                    atoms_per_frame = atom_count
                    
                # 读取注释行(可能包含盒子信息)
                line_index += 1
                comment = lines[line_index].strip()
                box_match = re.search(r'Box:\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', comment)
                if box_match:
                    box_size = np.array([
                        float(box_match.group(1)),
                        float(box_match.group(2)),
                        float(box_match.group(3))
                    ])
                
                # 读取原子坐标
                frame_elements = []
                frame_coords = []
                
                line_index += 1
                for i in range(atom_count):
                    if line_index >= len(lines):
                        break
                        
                    line = lines[line_index].strip()
                    parts = line.split()
                    
                    if len(parts) >= 4:
                        element = parts[0]
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        
                        frame_elements.append(element)
                        frame_coords.append([x, y, z])
                    
                    line_index += 1
                
                frames.append(np.array(frame_coords))
                if not elements:  # 只保存第一帧的元素列表，假设所有帧元素相同
                    elements = frame_elements
        
        if not frames:
            self.logger.error(f"未能从文件 {trajectory_file} 读取任何帧")
            return None, None, None
        
        # 选择指定帧
        if frame_index < 0:
            frame_index = len(frames) + frame_index
        
        if frame_index < 0 or frame_index >= len(frames):
            self.logger.error(f"帧索引 {frame_index} 超出范围 [0, {len(frames)-1}]")
            return None, None, None
        
        return elements, frames[frame_index], box_size
    
    def identify_molecules(self, elements, coords, bond_threshold=2.0):
        """识别分子结构
        
        Args:
            elements: 元素列表
            coords: 坐标数组
            bond_threshold: 成键距离阈值
            
        Returns:
            分子索引列表的列表（每个分子的原子索引）
        """
        self.logger.info("基于原子距离识别分子")
        
        n_atoms = len(elements)
        # 创建连接矩阵
        connections = np.zeros((n_atoms, n_atoms), dtype=bool)
        
        # 计算所有原子对之间的距离
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < bond_threshold:
                    connections[i, j] = connections[j, i] = True
        
        # 使用深度优先搜索找到连通分量（即分子）
        visited = np.zeros(n_atoms, dtype=bool)
        molecules = []
        
        for i in range(n_atoms):
            if not visited[i]:
                molecule = []
                self._dfs(i, connections, visited, molecule)
                molecules.append(molecule)
        
        self.logger.info(f"找到 {len(molecules)} 个分子")
        return molecules
    
    def _dfs(self, node, connections, visited, component):
        """深度优先搜索用于找到连通分量（分子）
        
        Args:
            node: 当前节点
            connections: 连接矩阵
            visited: 访问标记数组
            component: 当前连通分量
        """
        visited[node] = True
        component.append(node)
        
        for neighbor in np.where(connections[node])[0]:
            if not visited[neighbor]:
                self._dfs(neighbor, connections, visited, component)
    
    def identify_ion_types(self, elements, molecules):
        """识别分子类型（阳离子、阴离子、溶剂）
        
        Args:
            elements: 元素列表
            molecules: 分子索引列表的列表
            
        Returns:
            (阳离子列表, 阴离子列表, 溶剂列表)，每个列表包含分子索引
        """
        self.logger.info("识别分子类型（阳离子、阴离子、溶剂）")
        
        # 定义阳离子和阴离子元素
        cation_elements = ['Li', 'Na', 'K', 'Zn', 'Mg', 'Ca', 'Rb', 'Cs']
        anion_elements = ['F', 'Cl', 'Br', 'I', 'O', 'N', 'S', 'P', 'B']
        
        cations = []
        anions = []
        solvents = []
        
        for mol_indices in molecules:
            # 检查分子是否是单原子离子
            if len(mol_indices) == 1:
                element = elements[mol_indices[0]]
                if element in cation_elements:
                    cations.append(mol_indices)
                elif element in anion_elements:
                    anions.append(mol_indices)
                else:
                    solvents.append(mol_indices)
            else:
                # 多原子分子，检查元素构成
                mol_elements = [elements[i] for i in mol_indices]
                
                # 简单启发式规则：含有N、O、S、P的多原子分子可能是阴离子
                if any(e in anion_elements for e in mol_elements) and len(mol_indices) < 20:
                    anions.append(mol_indices)
                else:
                    # 其他多原子分子视为溶剂
                    solvents.append(mol_indices)
        
        self.logger.info(f"识别出 {len(cations)} 个阳离子, {len(anions)} 个阴离子, {len(solvents)} 个溶剂分子")
        return cations, anions, solvents
    
    def unwrap_coordinates(self, coords, molecules, box_size):
        """修复周期性边界条件，使分子不跨越边界
        
        Args:
            coords: 坐标数组
            molecules: 分子索引列表的列表
            box_size: 盒子尺寸
            
        Returns:
            修复后的坐标数组
        """
        self.logger.info("修复周期性边界条件")
        
        fixed_coords = np.copy(coords)
        
        for mol_indices in molecules:
            if len(mol_indices) <= 1:
                continue  # 单原子分子不需要修复
                
            # 提取分子坐标
            mol_coords = coords[mol_indices]
            
            # 计算分子质心
            centroid = np.mean(mol_coords, axis=0)
            
            # 将所有原子移至靠近质心的周期性映像
            for idx, atom_idx in enumerate(mol_indices):
                for dim in range(3):
                    # 计算原子与质心在各个维度的距离
                    dist = fixed_coords[atom_idx, dim] - centroid[dim]
                    # 如果距离大于半个盒子，则调整位置
                    if dist > box_size[dim] / 2:
                        fixed_coords[atom_idx, dim] -= box_size[dim]
                    elif dist < -box_size[dim] / 2:
                        fixed_coords[atom_idx, dim] += box_size[dim]
        
        return fixed_coords
    
    def calculate_ion_solvent_distances(self, coords, cations, solvents):
        """计算所有阳离子和溶剂分子之间的距离
        
        Args:
            coords: 坐标数组
            cations: 阳离子列表
            solvents: 溶剂列表
            
        Returns:
            距离矩阵[阳离子索引][溶剂索引] = 距离
        """
        self.logger.info("计算阳离子-溶剂距离")
        
        distances = {}
        
        for i, cation_indices in enumerate(cations):
            distances[i] = {}
            
            # 计算阳离子位置（单原子离子或多原子离子的质心）
            if len(cation_indices) == 1:
                cation_pos = coords[cation_indices[0]]
            else:
                cation_pos = np.mean(coords[cation_indices], axis=0)
            
            # 计算与每个溶剂分子的距离
            for j, solvent_indices in enumerate(solvents):
                # 计算溶剂分子质心
                solvent_centroid = np.mean(coords[solvent_indices], axis=0)
                
                # 计算欧几里得距离
                dist = np.linalg.norm(cation_pos - solvent_centroid)
                distances[i][j] = dist
        
        return distances
    
    def extract_solvation_shell(self, elements, coords, cations, solvents, distances, 
                               cutoff_distance=5.0, max_solvents=20):
        """提取阳离子周围的溶剂化壳层
        
        Args:
            elements: 元素列表
            coords: 坐标数组
            cations: 阳离子列表
            solvents: 溶剂列表
            distances: 距离矩阵
            cutoff_distance: 截断距离
            max_solvents: 每个阳离子最多考虑的溶剂分子数
            
        Returns:
            溶剂化壳层信息的列表，每个元素包含(阳离子类型, 阳离子索引, 溶剂索引列表)
        """
        self.logger.info(f"提取溶剂化壳层，截断距离: {cutoff_distance}埃，最大溶剂数: {max_solvents}")
        
        solvation_shells = []
        
        for i, cation_indices in enumerate(cations):
            # 获取阳离子元素类型
            cation_element = elements[cation_indices[0]]
            
            # 按距离排序溶剂分子
            sorted_solvents = sorted(distances[i].items(), key=lambda x: x[1])
            
            # 只保留截断距离内的溶剂分子，最多max_solvents个
            close_solvents = [j for j, dist in sorted_solvents 
                             if dist <= cutoff_distance][:max_solvents]
            
            if close_solvents:
                solvation_shells.append((cation_element, i, close_solvents))
        
        self.logger.info(f"提取了 {len(solvation_shells)} 个溶剂化壳层")
        return solvation_shells
    
    def write_solvation_shell(self, output_path, elements, coords, cation_indices, solvent_indices_list, cation_element, cation_index):
        """将溶剂化壳层写入XYZ文件
        
        Args:
            output_path: 输出文件路径
            elements: 元素列表
            coords: 坐标数组
            cation_indices: 阳离子原子索引列表
            solvent_indices_list: 溶剂分子索引列表的列表
            cation_element: 阳离子元素类型
            cation_index: 阳离子索引
            
        Returns:
            是否成功
        """
        # 收集所有原子的索引
        all_indices = cation_indices.copy()
        for solvent_indices in solvent_indices_list:
            all_indices.extend(solvent_indices)
        
        # 构建输出数据
        selected_elements = [elements[i] for i in all_indices]
        selected_coords = coords[all_indices]
        
        # 计算阳离子位置（作为参考点）
        if len(cation_indices) == 1:
            reference_point = coords[cation_indices[0]]
        else:
            reference_point = np.mean(coords[cation_indices], axis=0)
        
        # 平移坐标，使阳离子位于中心
        centered_coords = selected_coords - reference_point
        
        # 创建输出目录
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
        
        # 写入XYZ文件
        try:
            with open(output_path, 'w') as f:
                f.write(f"{len(selected_elements)}\n")
                f.write(f"{cation_element}{cation_index+1}的溶剂化结构, 包含{len(solvent_indices_list)}个溶剂分子\n")
                
                for i, element in enumerate(selected_elements):
                    x, y, z = centered_coords[i]
                    f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")
                    
            self.logger.info(f"已保存溶剂化壳层到文件: {output_path}")
            return True
        except Exception as e:
            self.logger.error(f"保存溶剂化壳层文件失败: {str(e)}")
            return False
    
    def analyze_solvation_structure(self, trajectory_file, output_dir, frame_index=-1, cutoff_distance=5.0, max_solvents=20, bond_threshold=2.0):
        """分析溶剂化结构并输出结果
        
        Args:
            trajectory_file: 轨迹文件路径
            output_dir: 输出目录路径
            frame_index: 帧索引，-1表示最后一帧
            cutoff_distance: 截断距离，用于定义溶剂化壳层
            max_solvents: 每个阳离子最多考虑的溶剂分子数
            bond_threshold: 识别分子时的键长阈值
            
        Returns:
            (溶剂化壳层信息列表, 输出文件路径列表)
        """
        self.logger.info(f"开始分析文件 {trajectory_file} 的溶剂化结构")
        
        # 创建输出目录
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 读取轨迹数据
        elements, coords, box_size = self.read_xyz_trajectory(trajectory_file, frame_index)
        if elements is None:
            self.logger.error("无法读取轨迹数据")
            return [], []
        
        # 识别分子
        molecules = self.identify_molecules(elements, coords, bond_threshold)
        
        # 识别离子类型
        cations, anions, solvents = self.identify_ion_types(elements, molecules)
        
        # 修复周期性边界条件
        fixed_coords = self.unwrap_coordinates(coords, molecules, box_size)
        
        # 计算离子-溶剂距离
        distances = self.calculate_ion_solvent_distances(fixed_coords, cations, solvents)
        
        # 提取溶剂化壳层
        solvation_shells = self.extract_solvation_shell(
            elements, fixed_coords, cations, solvents, distances, 
            cutoff_distance, max_solvents
        )
        
        # 输出溶剂化壳层到文件
        output_files = []
        for cation_element, cation_idx, solvent_indices in solvation_shells:
            # 使用阳离子类型和索引命名文件
            output_file = output_dir / f"{cation_element}_{cation_idx+1}_solvation.xyz"
            
            # 获取该阳离子的溶剂分子
            solvent_molecules = [solvents[j] for j in solvent_indices]
            
            # 写入溶剂化壳层
            if self.write_solvation_shell(
                output_file, elements, fixed_coords, 
                cations[cation_idx], solvent_molecules,
                cation_element, cation_idx
            ):
                output_files.append(output_file)
        
        self.logger.info(f"溶剂化结构分析完成，输出了 {len(output_files)} 个文件")
        return solvation_shells, output_files
    
    def analyze_coordination_numbers(self, trajectory_file, output_dir, cutoff_distances=None, 
                                    frame_indices=None, bond_threshold=2.0):
        """分析配位数随时间/距离的变化
        
        Args:
            trajectory_file: 轨迹文件路径
            output_dir: 输出目录路径
            cutoff_distances: 截断距离列表，默认为[3.0, 4.0, 5.0]
            frame_indices: 要分析的帧索引列表，默认为所有帧
            bond_threshold: 识别分子时的键长阈值
            
        Returns:
            配位数数据和图表文件路径
        """
        self.logger.info(f"分析文件 {trajectory_file} 的配位数")
        
        if cutoff_distances is None:
            cutoff_distances = [3.0, 4.0, 5.0]
        
        # 创建输出目录
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 读取所有帧
        frames = []
        elements = None
        box_size = None
        
        with open(trajectory_file, 'r') as f:
            lines = f.readlines()
            
            line_index = 0
            while line_index < len(lines):
                if not lines[line_index].strip():
                    line_index += 1
                    continue
                
                # 读取原子数
                atom_count = int(lines[line_index].strip())
                    
                # 读取注释行(可能包含盒子信息)
                line_index += 1
                comment = lines[line_index].strip()
                box_match = re.search(r'Box:\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', comment)
                if box_match and box_size is None:
                    box_size = np.array([
                        float(box_match.group(1)),
                        float(box_match.group(2)),
                        float(box_match.group(3))
                    ])
                
                # 读取原子坐标
                frame_elements = []
                frame_coords = []
                
                line_index += 1
                for i in range(atom_count):
                    if line_index >= len(lines):
                        break
                        
                    line = lines[line_index].strip()
                    parts = line.split()
                    
                    if len(parts) >= 4:
                        element = parts[0]
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        
                        frame_elements.append(element)
                        frame_coords.append([x, y, z])
                    
                    line_index += 1
                
                frames.append(np.array(frame_coords))
                if elements is None:  # 只保存第一帧的元素列表，假设所有帧元素相同
                    elements = frame_elements
        
        if not frames:
            self.logger.error(f"未能从文件 {trajectory_file} 读取任何帧")
            return None, []
        
        # 选择要分析的帧
        if frame_indices is None:
            frame_indices = list(range(len(frames)))
        
        # 分析每一帧
        coordination_data = {cutoff: [] for cutoff in cutoff_distances}
        
        for frame_idx in frame_indices:
            if frame_idx < 0 or frame_idx >= len(frames):
                continue
                
            coords = frames[frame_idx]
            
            # 分析该帧
            molecules = self.identify_molecules(elements, coords, bond_threshold)
            cations, anions, solvents = self.identify_ion_types(elements, molecules)
            fixed_coords = self.unwrap_coordinates(coords, molecules, box_size)
            distances = self.calculate_ion_solvent_distances(fixed_coords, cations, solvents)
            
            # 计算不同截断距离下的配位数
            for cutoff in cutoff_distances:
                coordinations = []
                for cation_idx in distances:
                    # 计算该阳离子的配位数
                    coordination = sum(1 for dist in distances[cation_idx].values() if dist <= cutoff)
                    coordinations.append(coordination)
                
                # 计算平均配位数
                if coordinations:
                    avg_coordination = sum(coordinations) / len(coordinations)
                    coordination_data[cutoff].append((frame_idx, avg_coordination))
        
        # 绘制配位数随时间的变化图
        plt.figure(figsize=(10, 6))
        for cutoff in cutoff_distances:
            if coordination_data[cutoff]:
                x, y = zip(*coordination_data[cutoff])
                plt.plot(x, y, label=f'Cutoff = {cutoff} Å')
        
        plt.xlabel('Frame Index')
        plt.ylabel('Average Coordination Number')
        plt.title('Coordination Number vs Frame Index')
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # 保存图表
        plot_file = output_dir / "coordination_number_vs_frame.png"
        plt.savefig(plot_file, dpi=300)
        plt.close()
        
        # 保存数据到CSV
        csv_file = output_dir / "coordination_data.csv"
        with open(csv_file, 'w') as f:
            f.write("Frame")
            for cutoff in cutoff_distances:
                f.write(f",Coordination (cutoff={cutoff})")
            f.write("\n")
            
            # 假设所有截断距离的数据点数相同
            if coordination_data[cutoff_distances[0]]:
                for i in range(len(coordination_data[cutoff_distances[0]])):
                    frame_idx = coordination_data[cutoff_distances[0]][i][0]
                    f.write(f"{frame_idx}")
                    
                    for cutoff in cutoff_distances:
                        if i < len(coordination_data[cutoff]):
                            coord = coordination_data[cutoff][i][1]
                            f.write(f",{coord:.4f}")
                        else:
                            f.write(",")
                    
                    f.write("\n")
        
        self.logger.info(f"配位数分析完成，输出数据到: {csv_file}")
        return coordination_data, [plot_file, csv_file]
    
    def process_trajectory(self, trajectory_file, output_dir, unwrap=True, bond_threshold=2.0):
        """处理轨迹文件，修复周期性边界条件并生成新的轨迹
        
        Args:
            trajectory_file: 轨迹文件路径
            output_dir: 输出目录路径
            unwrap: 是否展开分子，防止跨越周期性边界
            bond_threshold: 识别分子时的键长阈值
            
        Returns:
            (是否成功, 输出文件路径)
        """
        self.logger.info(f"处理轨迹文件 {trajectory_file}")
        
        # 使用trajectory_processor处理轨迹
        trajectory_path = Path(trajectory_file)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = output_dir / f"{trajectory_path.stem}_fixed.xyz"
        
        success = self.trajectory_processor.fix_periodic_molecules(
            str(trajectory_path),
            str(output_file),
            bond_threshold=bond_threshold,
            unwrap=unwrap
        )
        
        if success:
            self.logger.info(f"轨迹处理成功，输出到: {output_file}")
        else:
            self.logger.error(f"轨迹处理失败")
        
        return success, output_file
    
    def analyze_ion_pairs(self, coords, cations, anions, solvents, 
                         cip_distance=4.0, ssip_distance=7.0):
        """分析离子对类型：CIP（接触离子对）、SSIP（溶剂分离离子对）和AGG（聚集体）
        
        Args:
            coords: 坐标数组
            cations: 阳离子列表
            anions: 阴离子列表
            solvents: 溶剂列表
            cip_distance: CIP的最大距离阈值（埃）
            ssip_distance: SSIP的最大距离阈值（埃）
            
        Returns:
            包含离子对分析结果的字典
        """
        self.logger.info(f"分析离子对类型，CIP距离阈值: {cip_distance}埃，SSIP距离阈值: {ssip_distance}埃")
        
        # 初始化结果
        results = {
            'cip_count': 0,  # 接触离子对数量
            'ssip_count': 0, # 溶剂分离离子对数量
            'agg_count': 0,  # 聚集体数量
            'free_cation_count': 0, # 自由阳离子数量
            'cip_pairs': [],  # CIP离子对列表
            'ssip_pairs': [], # SSIP离子对列表
            'agg_groups': [],  # 聚集体列表
            'free_cations': [] # 自由阳离子列表
        }
        
        # 如果没有阳离子或阴离子，直接返回
        if not cations or not anions:
            self.logger.warning("没有找到阳离子或阴离子，无法分析离子对")
            return results
        
        # 计算所有阳离子和阴离子之间的距离
        cation_anion_distances = {}
        for i, cation_indices in enumerate(cations):
            # 计算阳离子位置（单原子离子或多原子离子的质心）
            if len(cation_indices) == 1:
                cation_pos = coords[cation_indices[0]]
            else:
                cation_pos = np.mean(coords[cation_indices], axis=0)
            
            cation_anion_distances[i] = {}
            
            for j, anion_indices in enumerate(anions):
                # 计算阴离子位置
                if len(anion_indices) == 1:
                    anion_pos = coords[anion_indices[0]]
                else:
                    anion_pos = np.mean(coords[anion_indices], axis=0)
                
                # 计算欧几里得距离
                dist = np.linalg.norm(cation_pos - anion_pos)
                cation_anion_distances[i][j] = dist
        
        # 构建离子连接图（用于识别聚集体）
        ion_connections = {}
        for i in range(len(cations)):
            ion_connections[i] = []
            
        # 分析每个阳离子
        for i, cation_indices in enumerate(cations):
            # 按距离排序阴离子
            sorted_anions = sorted(cation_anion_distances[i].items(), key=lambda x: x[1])
            
            # 检查最近的阴离子
            if sorted_anions:
                nearest_anion_idx, nearest_distance = sorted_anions[0]
                
                # 判断离子对类型
                if nearest_distance <= cip_distance:
                    # 接触离子对 (CIP)
                    results['cip_count'] += 1
                    results['cip_pairs'].append((i, nearest_anion_idx))
                    ion_connections[i].append(nearest_anion_idx)
                elif nearest_distance <= ssip_distance:
                    # 溶剂分离离子对 (SSIP)
                    results['ssip_count'] += 1
                    results['ssip_pairs'].append((i, nearest_anion_idx))
                    ion_connections[i].append(nearest_anion_idx)
                else:
                    # 自由阳离子
                    results['free_cation_count'] += 1
                    results['free_cations'].append(i)
            else:
                # 没有阴离子
                results['free_cation_count'] += 1
                results['free_cations'].append(i)
        
        # 识别聚集体（连接的离子对网络）
        visited = set()
        
        for i in range(len(cations)):
            if i in visited:
                continue
                
            if not ion_connections[i]:
                continue  # 跳过自由离子
                
            # 使用BFS找到连接的离子对网络
            agg_group = []
            queue = [i]
            visited.add(i)
            
            while queue:
                current = queue.pop(0)
                agg_group.append(current)
                
                for neighbor in ion_connections[current]:
                    for cation_idx, pairs in ion_connections.items():
                        if cation_idx not in visited and neighbor in pairs:
                            queue.append(cation_idx)
                            visited.add(cation_idx)
            
            # 如果网络包含多个离子对，则为聚集体
            if len(agg_group) > 1:
                results['agg_count'] += 1
                results['agg_groups'].append(agg_group)
                
                # 从CIP和SSIP计数中减去已经计入聚集体的离子对
                for cation_idx in agg_group:
                    if any((cation_idx, anion_idx) in results['cip_pairs'] for anion_idx in range(len(anions))):
                        results['cip_count'] -= 1
                    elif any((cation_idx, anion_idx) in results['ssip_pairs'] for anion_idx in range(len(anions))):
                        results['ssip_count'] -= 1
        
        self.logger.info(f"离子对分析结果: CIP={results['cip_count']}, SSIP={results['ssip_count']}, AGG={results['agg_count']}, 自由阳离子={results['free_cation_count']}")
        return results
        
    def analyze_solvent_distribution(self, trajectory_file, atom_counts=None, 
                                    cation_atom_types=None, anion_atom_types=None, 
                                    solvent_atom_types=None, cutoff_distance=5.0,
                                    cip_distance=4.0, ssip_distance=7.0):
        """分析溶剂分布
        
        Args:
            trajectory_file: LAMMPS轨迹文件路径
            atom_counts: 分子原子数量映射
            cation_atom_types: 阳离子中代表性原子的类型列表
            anion_atom_types: 阴离子中代表性原子的类型列表
            solvent_atom_types: 溶剂中代表性原子的类型列表
            cutoff_distance: 溶剂化层的截止距离（埃）
            cip_distance: CIP的最大距离阈值（埃）
            ssip_distance: SSIP的最大距离阈值（埃）
            
        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"分析溶剂分布: {trajectory_file}")
        
        # 读取轨迹文件
        elements, coords, box_size = self.read_xyz_trajectory(trajectory_file)
        if elements is None:
            self.logger.error("无法读取轨迹数据")
            return None
        
        # 识别分子
        molecules = self.identify_molecules(elements, coords)
        
        # 识别离子类型
        cations, anions, solvents = self.identify_ion_types(elements, molecules)
        
        # 修复周期性边界条件
        fixed_coords = self.unwrap_coordinates(coords, molecules, box_size)
        
        # 计算离子-溶剂距离
        distances = self.calculate_ion_solvent_distances(fixed_coords, cations, solvents)
        
        # 分析离子对类型
        ion_pair_results = self.analyze_ion_pairs(fixed_coords, cations, anions, solvents, 
                                                 cip_distance, ssip_distance)
        
        # 计算每个阳离子周围的溶剂分布
        cation_solvation = {}
        
        for i, cation_indices in enumerate(cations):
            # 按距离排序溶剂分子
            sorted_solvents = sorted(distances[i].items(), key=lambda x: x[1])
            
            # 只考虑截断距离内的溶剂分子
            close_solvents = [(j, dist) for j, dist in sorted_solvents if dist <= cutoff_distance]
            
            # 统计溶剂分子数量
            cation_solvation[i] = len(close_solvents)
        
        # 计算平均溶剂化数
        avg_solvation = sum(cation_solvation.values()) / len(cation_solvation) if cation_solvation else 0
        
        # 统计溶剂化数分布
        solvation_distribution = {}
        for count in cation_solvation.values():
            solvation_distribution[count] = solvation_distribution.get(count, 0) + 1
        
        # 整合结果
        results = {
            'cation_count': len(cations),
            'anion_count': len(anions),
            'solvent_count': len(solvents),
            'avg_solvation': avg_solvation,
            'solvation_distribution': solvation_distribution,
            'cation_solvation': cation_solvation,
            'avg_cation_solvation': {},
            'frames': 1,
            'ion_pair_analysis': ion_pair_results
        }
        
        # 将离子对分析结果直接添加到结果字典中
        results['cip_count'] = ion_pair_results['cip_count']
        results['ssip_count'] = ion_pair_results['ssip_count']
        results['agg_count'] = ion_pair_results['agg_count']
        results['free_cation_count'] = ion_pair_results['free_cation_count']
        
        # 计算百分比
        total_cations = len(cations)
        if total_cations > 0:
            results['cip_percentage'] = (ion_pair_results['cip_count'] / total_cations) * 100
            results['ssip_percentage'] = (ion_pair_results['ssip_count'] / total_cations) * 100
            results['agg_percentage'] = (ion_pair_results['agg_count'] / total_cations) * 100
            results['free_cation_percentage'] = (ion_pair_results['free_cation_count'] / total_cations) * 100
        else:
            results['cip_percentage'] = 0
            results['ssip_percentage'] = 0
            results['agg_percentage'] = 0
            results['free_cation_percentage'] = 0
        
        self.logger.info(f"溶剂分布分析完成: 平均溶剂化数={avg_solvation:.2f}")
        return results

class TemperatureSeriesAnalyzer:
    def __init__(self, base_dir, formulation, temperatures):
        self.base_dir = base_dir
        self.formulation = formulation
        self.temperatures = temperatures
        self.results = {}
        
    def analyze_diffusion_temperature_dependence(self):
        """分析扩散系数的温度依赖性，计算活化能"""
        diffusion_data = []
        for temp in self.temperatures:
            sample_name = f"{self.formulation}_{temp}K"
            result = self.analyzer.run_msd_analysis(sample_name)
            diffusion_data.append((temp, result['diffusion_coeffs']['D_total']))
            
        # 计算活化能
        activation_energy = self.calculate_activation_energy(diffusion_data)
        # 绘制Arrhenius图
        arrhenius_plot = self.plot_arrhenius(diffusion_data, activation_energy)
        
        return {
            'diffusion_data': diffusion_data,
            'activation_energy': activation_energy,
            'arrhenius_plot': arrhenius_plot
        } 