"""
轨迹处理器模块，用于处理和修复分子动力学模拟轨迹
"""
import os
import numpy as np
from pathlib import Path
import re
from ..utils.logger import Logger

class TrajectoryProcessor:
    """轨迹处理器类，用于处理和修复分子动力学模拟轨迹"""
    
    def __init__(self):
        """初始化轨迹处理器"""
        self.logger = Logger().get_logger()
    
    def load_xyz_file(self, xyz_file_path):
        """加载XYZ文件
        
        Args:
            xyz_file_path: XYZ文件路径
            
        Returns:
            (原子数量, 注释, 原子类型列表, 坐标数组)
        """
        self.logger.info(f"加载XYZ文件: {xyz_file_path}")
        
        with open(xyz_file_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            self.logger.error(f"XYZ文件格式错误: {xyz_file_path}")
            return None, None, None, None
        
        # 读取原子数量
        n_atoms = int(lines[0].strip())
        # 读取注释行
        comment = lines[1].strip()
        
        # 解析盒子尺寸（如果在注释行中）
        box_size = None
        box_match = re.search(r'Box:\s*(\d+\.?\d*)\s*(\d+\.?\d*)\s*(\d+\.?\d*)', comment)
        if box_match:
            box_size = np.array([float(box_match.group(1)), 
                                 float(box_match.group(2)), 
                                 float(box_match.group(3))])
            self.logger.info(f"从XYZ文件中提取盒子尺寸: {box_size}")
        
        # 读取原子类型和坐标
        atom_types = []
        coords = np.zeros((n_atoms, 3))
        
        for i in range(n_atoms):
            if i + 2 >= len(lines):
                self.logger.error(f"XYZ文件行数不足: {xyz_file_path}")
                return None, None, None, None
                
            parts = lines[i + 2].strip().split()
            if len(parts) < 4:
                self.logger.error(f"XYZ文件格式错误，行 {i+3}: {lines[i+2]}")
                return None, None, None, None
                
            atom_types.append(parts[0])
            coords[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
        
        return n_atoms, comment, atom_types, coords, box_size
    
    def write_xyz_file(self, xyz_file_path, n_atoms, comment, atom_types, coords):
        """写入XYZ文件
        
        Args:
            xyz_file_path: 输出XYZ文件路径
            n_atoms: 原子数量
            comment: 注释行
            atom_types: 原子类型列表
            coords: 坐标数组
            
        Returns:
            是否成功
        """
        self.logger.info(f"写入XYZ文件: {xyz_file_path}")
        
        try:
            with open(xyz_file_path, 'w') as f:
                f.write(f"{n_atoms}\n")
                f.write(f"{comment}\n")
                
                for i in range(n_atoms):
                    f.write(f"{atom_types[i]:2s} {coords[i, 0]:12.6f} {coords[i, 1]:12.6f} {coords[i, 2]:12.6f}\n")
                    
            return True
        except Exception as e:
            self.logger.error(f"写入XYZ文件失败: {str(e)}")
            return False
    
    def get_molecule_indices(self, atom_types, coords, bond_threshold=2.0):
        """根据原子间距离确定分子的索引分组
        
        Args:
            atom_types: 原子类型列表
            coords: 坐标数组
            bond_threshold: 成键距离阈值（埃）
            
        Returns:
            分子索引列表的列表，每个子列表包含一个分子的所有原子索引
        """
        self.logger.info(f"根据原子间距离分析分子结构")
        
        n_atoms = len(atom_types)
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
    
    def wrap_coordinates(self, coords, box_size):
        """将坐标包装到模拟盒子内
        
        Args:
            coords: 原子坐标数组
            box_size: 盒子尺寸数组
            
        Returns:
            包装后的坐标数组
        """
        self.logger.info(f"将坐标包装到盒子内")
        
        wrapped_coords = np.copy(coords)
        wrapped_coords = wrapped_coords % box_size
        
        return wrapped_coords
    
    def unwrap_molecule(self, mol_coords, box_size):
        """展开单个分子的坐标，使其不跨越周期性边界
        
        Args:
            mol_coords: 分子坐标数组
            box_size: 盒子尺寸数组
            
        Returns:
            展开后的分子坐标数组
        """
        # 找到分子的质心
        centroid = np.mean(mol_coords, axis=0)
        
        # 将所有原子移至靠近质心的周期性映像
        unwrapped_coords = np.copy(mol_coords)
        for i in range(len(mol_coords)):
            for dim in range(3):
                # 计算原子与质心在各个维度的距离
                dist = unwrapped_coords[i, dim] - centroid[dim]
                # 如果距离大于半个盒子，则调整位置
                if dist > box_size[dim] / 2:
                    unwrapped_coords[i, dim] -= box_size[dim]
                elif dist < -box_size[dim] / 2:
                    unwrapped_coords[i, dim] += box_size[dim]
        
        return unwrapped_coords
    
    def fix_periodic_molecules(self, xyz_file_path, output_file_path=None, bond_threshold=2.0, unwrap=True):
        """修复由于周期性边界条件导致的分子断裂问题
        
        Args:
            xyz_file_path: 输入XYZ文件路径
            output_file_path: 输出XYZ文件路径，如果为None则在原文件名基础上添加后缀
            bond_threshold: 原子成键的距离阈值
            unwrap: 是否展开分子（True）或仅将所有分子移至盒子内（False）
            
        Returns:
            是否成功
        """
        self.logger.info(f"处理XYZ文件中的周期性边界问题: {xyz_file_path}")
        
        # 如果没有指定输出文件，则在原文件名基础上添加后缀
        if output_file_path is None:
            input_path = Path(xyz_file_path)
            output_file_path = input_path.parent / f"{input_path.stem}_fixed{input_path.suffix}"
        
        # 加载XYZ文件
        n_atoms, comment, atom_types, coords, box_size = self.load_xyz_file(xyz_file_path)
        
        if n_atoms is None or box_size is None:
            self.logger.error("加载XYZ文件失败或无法确定盒子尺寸")
            return False
        
        # 初始包装所有坐标到盒子内
        coords = self.wrap_coordinates(coords, box_size)
        
        if unwrap:
            # 获取分子索引分组
            molecules = self.get_molecule_indices(atom_types, coords, bond_threshold)
            
            # 处理每个分子，避免跨越周期性边界
            fixed_coords = np.copy(coords)
            for molecule in molecules:
                mol_coords = coords[molecule]
                # 展开分子坐标
                unwrapped_mol_coords = self.unwrap_molecule(mol_coords, box_size)
                # 更新坐标
                fixed_coords[molecule] = unwrapped_mol_coords
            
            self.logger.info(f"已修复 {len(molecules)} 个分子的周期性问题")
        else:
            # 仅包装所有坐标到盒子内
            fixed_coords = coords
            self.logger.info("已将所有坐标包装到盒子内")
        
        # 写入修复后的XYZ文件
        success = self.write_xyz_file(output_file_path, n_atoms, comment, atom_types, fixed_coords)
        
        if success:
            self.logger.info(f"已成功保存修复后的XYZ文件: {output_file_path}")
        
        return success
    
    def process_trajectory_directory(self, directory_path, file_pattern="*.xyz", unwrap=True):
        """处理目录中所有轨迹文件
        
        Args:
            directory_path: 目录路径
            file_pattern: 文件匹配模式
            unwrap: 是否展开分子
            
        Returns:
            处理成功的文件数量
        """
        self.logger.info(f"批量处理目录中的轨迹文件: {directory_path}")
        
        directory = Path(directory_path)
        if not directory.exists() or not directory.is_dir():
            self.logger.error(f"目录不存在: {directory_path}")
            return 0
        
        # 创建输出目录
        output_dir = directory / "fixed_trajectories"
        output_dir.mkdir(exist_ok=True)
        
        # 获取所有匹配的文件
        xyz_files = list(directory.glob(file_pattern))
        success_count = 0
        
        for xyz_file in xyz_files:
            try:
                output_file = output_dir / f"{xyz_file.stem}_fixed{xyz_file.suffix}"
                if self.fix_periodic_molecules(str(xyz_file), str(output_file), unwrap=unwrap):
                    success_count += 1
            except Exception as e:
                self.logger.error(f"处理文件 {xyz_file} 时出错: {str(e)}")
        
        self.logger.info(f"成功处理 {success_count}/{len(xyz_files)} 个文件")
        return success_count 