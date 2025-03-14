"""
高斯计算分析模块

用于处理和分析高斯量子化学计算的结果
"""

import os
import json
import logging
from typing import Dict, List, Tuple, Optional, Union, Any
import numpy as np

class GaussianAnalyzer:
    """
    高斯计算分析器
    
    提供处理和分析高斯量子化学计算结果的功能
    """
    def __init__(self, working_dir: str, energy_cache_path: Optional[str] = None, logger=None):
        """
        初始化高斯分析器

        Args:
            working_dir: 工作目录路径
            energy_cache_path: 能量缓存文件路径，如果为None则使用默认路径
            logger: 日志记录器实例
        """
        self.working_dir = working_dir
        self.energy_cache_path = energy_cache_path or os.path.join(working_dir, 'energy_cache.json')
        self.logger = logger or logging.getLogger(__name__)
    
    def extract_gaussian_energy(self, file_path: str) -> Optional[float]:
        """
        从Gaussian输出文件中提取能量，如果能量已经在缓存中，则直接返回

        Args:
            file_path: 高斯输出文件路径

        Returns:
            提取的能量值，如果无法提取则返回None
        """
        self.logger.info(f"从文件中提取高斯能量: {file_path}")
        
        # 检查缓存中是否已经存在该文件的能量
        energy_cache = self.load_energy_cache()
        if file_path in energy_cache:
            self.logger.debug(f"从缓存中获取能量: {file_path} -> {energy_cache[file_path]}")
            return energy_cache[file_path]

        if not os.path.exists(file_path):
            self.logger.warning(f"高斯输出文件不存在: {file_path}")
            return None

        try:
            # 从Gaussian文件中提取能量的逻辑
            with open(file_path, 'r') as f:
                for line in f:
                    if "SCF Done" in line:
                        energy = float(line.split()[4])
                        self.logger.debug(f"从文件中提取能量: {energy}")
                        
                        # 更新缓存
                        self.update_energy_cache(file_path, energy)
                        return energy
            
            self.logger.warning(f"在文件中未找到能量信息: {file_path}")
            return None
            
        except Exception as e:
            self.logger.error(f"提取高斯能量时出错: {e}")
            return None

    def is_gaussian_finished(self, file_path: str) -> bool:
        """
        检查Gaussian的log文件是否包含计算完成的标志

        Args:
            file_path: 高斯输出文件路径

        Returns:
            如果文件包含计算完成的标志，则返回True
        """
        self.logger.info(f"检查高斯计算是否完成: {file_path}")
        
        if not os.path.exists(file_path):
            self.logger.warning(f"高斯输出文件不存在: {file_path}")
            return False
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if "Normal termination" in line:
                        self.logger.debug(f"高斯计算正常完成")
                        return True
            
            self.logger.warning(f"高斯计算未正常完成")
            return False
            
        except Exception as e:
            self.logger.error(f"检查高斯计算完成状态时出错: {e}")
            return False

    def generate_gaussian_input(self, 
                                molecule_name: str, 
                                atom_positions: List[Tuple[str, Tuple[float, float, float]]], 
                                charge: int, 
                                spin: int, 
                                output_path: str,
                                method: str = 'B3LYP',
                                basis_set: str = '6-31G(d)',
                                nproc: int = 4,
                                mem: str = '4GB',
                                additional_keywords: Optional[List[str]] = None) -> bool:
        """
        生成Gaussian输入文件

        Args:
            molecule_name: 分子名称
            atom_positions: 原子位置列表，每个元素为(原子名称, (x, y, z))
            charge: 分子电荷
            spin: 自旋多重度
            output_path: 输出文件路径
            method: 计算方法
            basis_set: 基组
            nproc: 处理器数量
            mem: 内存大小
            additional_keywords: 额外的关键字参数

        Returns:
            如果文件生成成功，则返回True
        """
        self.logger.info(f"生成高斯输入文件: {output_path}")
        
        try:
            keywords = ['Opt', 'Freq']
            if additional_keywords:
                keywords.extend(additional_keywords)
            
            with open(output_path, 'w') as f:
                f.write(f"%nprocshared={nproc}\n")
                f.write(f"%mem={mem}\n")
                f.write(f"#P {method}/{basis_set} {' '.join(keywords)}\n\n")
                f.write(f"{molecule_name} Structure\n\n")
                f.write(f"{charge} {spin}\n")
                for atom_name, pos in atom_positions:
                    f.write(f"{atom_name} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
                f.write("\n")
            
            self.logger.info(f"高斯输入文件已生成: {output_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成高斯输入文件时出错: {e}")
            return False

    def load_energy_cache(self) -> Dict[str, float]:
        """
        从指定路径加载能量缓存

        Returns:
            能量缓存字典，如果缓存文件不存在，返回空字典
        """
        if os.path.exists(self.energy_cache_path):
            try:
                with open(self.energy_cache_path, 'r') as cache_file:
                    cache = json.load(cache_file)
                self.logger.debug(f"已加载能量缓存，包含 {len(cache)} 个条目")
                return cache
            except Exception as e:
                self.logger.error(f"加载能量缓存时出错: {e}")
                return {}
        else:
            self.logger.debug(f"能量缓存文件不存在，返回空字典")
            return {}

    def update_energy_cache(self, file_path: str, energy: float) -> None:
        """
        更新能量缓存文件，保存新计算的能量

        Args:
            file_path: 高斯输出文件路径
            energy: 能量值
        """
        self.logger.info(f"更新能量缓存: {file_path} -> {energy}")
        
        try:
            # 加载现有缓存
            energy_cache = self.load_energy_cache()
            
            # 更新缓存
            energy_cache[file_path] = energy
            
            # 保存缓存
            with open(self.energy_cache_path, 'w') as cache_file:
                json.dump(energy_cache, cache_file, indent=2)
            
            self.logger.debug(f"能量缓存已更新，现包含 {len(energy_cache)} 个条目")
            
        except Exception as e:
            self.logger.error(f"更新能量缓存时出错: {e}")
    
    def analyze_all_gaussian_outputs(self, directory: str) -> Dict[str, Dict[str, Any]]:
        """
        分析指定目录下的所有高斯输出文件

        Args:
            directory: 包含高斯输出文件的目录路径

        Returns:
            包含每个文件分析结果的字典
        """
        self.logger.info(f"分析目录中的所有高斯输出: {directory}")
        
        results = {}
        
        try:
            # 遍历目录中的所有.log文件
            for filename in os.listdir(directory):
                if filename.endswith('.log'):
                    file_path = os.path.join(directory, filename)
                    
                    # 检查计算是否完成
                    is_finished = self.is_gaussian_finished(file_path)
                    
                    # 提取能量（如果计算完成）
                    energy = None
                    if is_finished:
                        energy = self.extract_gaussian_energy(file_path)
                    
                    # 记录结果
                    results[filename] = {
                        'path': file_path,
                        'finished': is_finished,
                        'energy': energy,
                        'basename': os.path.splitext(filename)[0]
                    }
            
            self.logger.info(f"高斯输出分析完成，分析了 {len(results)} 个文件")
            return results
            
        except Exception as e:
            self.logger.error(f"分析高斯输出时出错: {e}")
            return results
    
    def extract_vibrational_frequencies(self, file_path: str) -> List[float]:
        """
        从高斯输出文件中提取振动频率

        Args:
            file_path: 高斯输出文件路径

        Returns:
            振动频率列表
        """
        self.logger.info(f"从文件中提取振动频率: {file_path}")
        
        frequencies = []
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                
                for i, line in enumerate(lines):
                    if "Frequencies" in line:
                        # 高斯输出中的频率行格式类似于：
                        # Frequencies --   123.4567   234.5678   345.6789
                        parts = line.split()
                        for part in parts[2:]:  # 跳过 "Frequencies" 和 "--"
                            try:
                                freq = float(part)
                                frequencies.append(freq)
                            except ValueError:
                                continue
            
            self.logger.debug(f"提取到 {len(frequencies)} 个振动频率")
            return frequencies
            
        except Exception as e:
            self.logger.error(f"提取振动频率时出错: {e}")
            return []
    
    def calculate_thermochemistry(self, file_path: str) -> Dict[str, float]:
        """
        从高斯输出文件中提取热力学数据

        Args:
            file_path: 高斯输出文件路径

        Returns:
            包含热力学数据的字典
        """
        self.logger.info(f"从文件中提取热力学数据: {file_path}")
        
        thermo_data = {}
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                
                for i, line in enumerate(lines):
                    # 零点能
                    if "Zero-point correction" in line:
                        thermo_data['zero_point_energy'] = float(line.split("=")[1].split()[0])
                    
                    # 热校正到能量
                    elif "Thermal correction to Energy" in line:
                        thermo_data['thermal_energy'] = float(line.split("=")[1].strip())
                    
                    # 热校正到焓
                    elif "Thermal correction to Enthalpy" in line:
                        thermo_data['thermal_enthalpy'] = float(line.split("=")[1].strip())
                    
                    # 热校正到吉布斯自由能
                    elif "Thermal correction to Gibbs Free Energy" in line:
                        thermo_data['thermal_free_energy'] = float(line.split("=")[1].strip())
                    
                    # 总能量（计算水平下的能量）
                    elif "Sum of electronic and zero-point Energies" in line:
                        thermo_data['total_energy_with_zpe'] = float(line.split("=")[1].strip())
                    
                    # 总能量+热能
                    elif "Sum of electronic and thermal Energies" in line:
                        thermo_data['total_thermal_energy'] = float(line.split("=")[1].strip())
                    
                    # 总能量+焓
                    elif "Sum of electronic and thermal Enthalpies" in line:
                        thermo_data['total_thermal_enthalpy'] = float(line.split("=")[1].strip())
                    
                    # 总能量+吉布斯自由能
                    elif "Sum of electronic and thermal Free Energies" in line:
                        thermo_data['total_free_energy'] = float(line.split("=")[1].strip())
            
            self.logger.debug(f"提取到的热力学数据: {thermo_data}")
            return thermo_data
            
        except Exception as e:
            self.logger.error(f"提取热力学数据时出错: {e}")
            return {} 