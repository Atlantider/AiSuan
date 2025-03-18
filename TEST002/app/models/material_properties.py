import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.io.xyz import XYZ
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints, Potcar
import os
import json
import traceback

class MaterialPropertiesCalculator:
    """材料属性计算器"""
    
    def __init__(self):
        """初始化计算器"""
        pass
    
    def calculate(self, data):
        """
        计算材料属性
        
        参数:
        data (dict): 包含计算参数的字典
            - structure_file_path: 结构文件路径
            - properties: 要计算的属性列表
            
        返回:
        dict: 计算结果
        """
        try:
            # 解析参数
            structure_path = data.get('structure_file_path')
            properties = data.get('properties', [])
            
            print(f"开始计算材料属性，文件路径: {structure_path}")
            print(f"要计算的属性: {properties}")
            
            # 检查文件是否存在
            if not os.path.exists(structure_path):
                raise FileNotFoundError(f"结构文件不存在: {structure_path}")
            
            # 加载结构
            try:
                structure = Structure.from_file(structure_path)
                print(f"成功加载结构: {structure}")
            except Exception as e:
                print(f"加载结构文件时出错: {str(e)}")
                raise ValueError(f"无法解析结构文件: {str(e)}")
            
            # 初始化结果字典
            result = {}
            
            # 计算基本属性
            result['formula'] = structure.composition.reduced_formula
            print(f"化学式: {result['formula']}")
            
            # 计算密度
            if 'density' in properties:
                result['density'] = structure.density
                print(f"密度: {result['density']} g/cm³")
            
            # 计算晶格参数
            if 'structure' in properties:
                lattice = structure.lattice
                result['lattice_parameters'] = {
                    'a': lattice.a,
                    'b': lattice.b,
                    'c': lattice.c,
                    'alpha': lattice.alpha,
                    'beta': lattice.beta,
                    'gamma': lattice.gamma
                }
                print(f"晶格参数: a={lattice.a}, b={lattice.b}, c={lattice.c}")
                
                # 计算空间群
                try:
                    spg_analyzer = SpacegroupAnalyzer(structure)
                    result['space_group'] = spg_analyzer.get_space_group_symbol()
                    print(f"空间群: {result['space_group']}")
                except Exception as e:
                    print(f"计算空间群时出错: {str(e)}")
                    result['space_group'] = "无法确定"
            
            # 无论是否选择了structure属性，都生成结构数据用于可视化
            try:
                # 使用自定义方法生成XYZ格式数据
                xyz_data = self._generate_xyz_data(structure)
                result['structure_data'] = xyz_data
                print(f"生成XYZ数据，长度: {len(result['structure_data'])}")
                print(f"XYZ数据前100个字符: {result['structure_data'][:100]}")
            except Exception as e:
                print(f"生成XYZ数据时出错: {str(e)}")
                print(traceback.format_exc())
                # 生成简单的XYZ数据作为备用
                result['structure_data'] = self._generate_fallback_xyz(structure)
            
            # 计算能带结构（模拟数据，实际应用中应使用DFT计算）
            if 'bandstructure' in properties:
                print("生成能带结构数据")
                # 这里使用模拟数据，实际应用中应从DFT计算获取
                result['bandstructure_data'] = self._generate_mock_bandstructure(structure)
            
            # 计算态密度（模拟数据，实际应用中应使用DFT计算）
            if 'dos' in properties:
                print("生成态密度数据")
                # 这里使用模拟数据，实际应用中应从DFT计算获取
                result['dos_data'] = self._generate_mock_dos(structure)
            
            # 新增: 计算体积
            if 'volume' in properties:
                print("计算晶胞体积")
                result['volume'] = self._calculate_volume(structure)
            
            # 新增: 计算最短原子距离
            if 'min_distance' in properties:
                print("计算最短原子距离")
                result['min_distance'] = self._calculate_min_distance(structure)
            
            # 新增: 计算表面能（估算值）
            if 'surface_energy' in properties:
                print("估算表面能")
                result['surface_energy'] = self._estimate_surface_energies(structure)
            
            # 新增: 生成VASP输入文件
            if 'vasp_input' in properties:
                print("生成VASP输入文件")
                result['vasp_input'] = self._generate_vasp_input(structure)
            
            print("材料属性计算完成")
            return result
            
        except Exception as e:
            print(f"计算材料属性时出错: {str(e)}")
            print(traceback.format_exc())
            raise Exception(f"计算材料属性时出错: {str(e)}")
    
    def _generate_xyz_data(self, structure):
        """
        生成XYZ格式的结构数据
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        str: XYZ格式的结构数据
        """
        # 获取原子数量
        num_atoms = len(structure)
        
        # 构建XYZ格式字符串
        xyz_lines = [str(num_atoms), structure.composition.reduced_formula]
        
        # 添加每个原子的信息
        for site in structure:
            # 获取元素符号
            symbol = site.species_string
            # 获取坐标
            x, y, z = site.coords
            # 添加到XYZ数据中
            xyz_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
        
        # 合并为一个字符串
        xyz_data = "\n".join(xyz_lines)
        return xyz_data
    
    def _generate_fallback_xyz(self, structure):
        """
        生成简单的XYZ备用数据
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        str: 简单的XYZ格式数据
        """
        # 创建一个简单的立方体结构
        formula = structure.composition.reduced_formula
        num_atoms = min(8, len(structure))  # 最多取8个原子
        
        xyz_lines = [str(num_atoms), f"Fallback structure for {formula}"]
        
        # 添加立方体顶点的原子
        symbols = [site.species_string for site in structure][:num_atoms]
        positions = [
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
            (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)
        ][:num_atoms]
        
        for i in range(num_atoms):
            symbol = symbols[i]
            x, y, z = positions[i]
            xyz_lines.append(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
        
        return "\n".join(xyz_lines)
    
    def _generate_mock_bandstructure(self, structure):
        """
        生成模拟的能带结构数据
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 能带结构数据
        """
        try:
            # 模拟k点路径
            kpoints = np.linspace(0, 1, 100)
            
            # 模拟能带数据（实际应从DFT计算获取）
            # 这里简单地生成几条正弦波形的能带
            num_bands = 5
            energies = []
            
            for i in range(num_bands):
                # 生成一条能带
                band = 2 * np.sin(2 * np.pi * kpoints) + i * 2 - 4
                energies.append(band.tolist())
            
            # 高对称点标记
            tick_vals = [0, 0.33, 0.66, 1]
            tick_labels = ['Γ', 'X', 'M', 'Γ']
            
            return {
                'kpoints': kpoints.tolist(),
                'energies': energies,
                'tick_vals': tick_vals,
                'tick_labels': tick_labels
            }
        except Exception as e:
            print(f"生成模拟能带数据时出错: {str(e)}")
            # 返回简化的数据
            return {
                'kpoints': [0, 0.5, 1],
                'energies': [[-1, 0, -1], [1, 2, 1]],
                'tick_vals': [0, 0.5, 1],
                'tick_labels': ['Γ', 'X', 'Γ']
            }
    
    def _generate_mock_dos(self, structure):
        """
        生成模拟的态密度数据
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 态密度数据
        """
        try:
            # 模拟能量范围
            energies = np.linspace(-10, 10, 200)
            
            # 模拟态密度（实际应从DFT计算获取）
            # 这里简单地生成高斯分布的态密度
            densities = np.exp(-(energies)**2/5) + 0.5 * np.exp(-(energies-5)**2/2)
            
            return {
                'energies': energies.tolist(),
                'densities': densities.tolist()
            }
        except Exception as e:
            print(f"生成模拟态密度数据时出错: {str(e)}")
            # 返回简化的数据
            return {
                'energies': [-5, 0, 5],
                'densities': [0.1, 1.0, 0.1]
            }
    
    # 新增方法: 计算晶胞体积
    def _calculate_volume(self, structure):
        """
        计算晶胞体积
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 体积相关数据
        """
        try:
            # 计算晶胞体积
            volume = structure.volume
            
            # 计算每个原子的体积
            volume_per_atom = volume / len(structure)
            
            # 计算每个公式单位的体积
            formula_units = structure.composition.get_reduced_formula_and_factor()[1]
            volume_per_formula = volume / formula_units
            
            return {
                'cell_volume': volume,  # 晶胞体积 (Å³)
                'volume_per_atom': volume_per_atom,  # 每原子体积 (Å³/atom)
                'volume_per_formula': volume_per_formula,  # 每公式单位体积 (Å³/formula)
                'formula_units': formula_units  # 晶胞中的公式单位数
            }
        except Exception as e:
            print(f"计算体积时出错: {str(e)}")
            return {
                'cell_volume': structure.volume,
                'error': str(e)
            }
    
    # 新增方法: 计算最短原子距离
    def _calculate_min_distance(self, structure):
        """
        计算结构中的最短原子距离和配位数
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 原子距离和配位数据
        """
        try:
            print("开始计算最短原子距离和配位数")
            
            # 计算最短原子距离
            min_distance = float('inf')
            for i, site1 in enumerate(structure):
                for j, site2 in enumerate(structure):
                    if i != j:  # 不比较同一个原子
                        distance = structure.get_distance(i, j)
                        min_distance = min(min_distance, distance)
            
            print(f"最短原子距离: {min_distance} Å")
            
            # 使用CrystalNN计算配位环境
            print("初始化CrystalNN分析器")
            nn_analyzer = CrystalNN()
            coordination_numbers = []
            coordination_details = []
            
            print(f"结构中共有 {len(structure)} 个原子")
            for i, site in enumerate(structure):
                try:
                    print(f"计算原子 {i} ({site.species_string}) 的配位环境")
                    # 获取近邻信息
                    nn_info = nn_analyzer.get_nn_info(structure, i)
                    
                    # 配位数
                    cn = len(nn_info)
                    coordination_numbers.append(cn)
                    print(f"原子 {i} 的配位数: {cn}")
                    
                    # 配位细节
                    neighbors = []
                    for neighbor in nn_info:
                        try:
                            # 适应可能的API变化
                            if isinstance(neighbor, dict):
                                # 旧版API返回字典
                                element = neighbor['site'].species_string if hasattr(neighbor['site'], 'species_string') else str(neighbor['site'].species)
                                distance = neighbor.get('distance')
                                weight = neighbor.get('weight', 1.0)
                            else:
                                # 新版API可能返回对象
                                element = neighbor.site.species_string if hasattr(neighbor.site, 'species_string') else str(neighbor.site.species)
                                distance = getattr(neighbor, 'distance', None)
                                weight = getattr(neighbor, 'weight', 1.0)
                            
                            neighbors.append({
                                'element': element,
                                'distance': distance,
                                'weight': weight
                            })
                            print(f"  - 近邻: {element}, 距离: {distance} Å, 权重: {weight}")
                        except Exception as e:
                            print(f"处理近邻信息时出错: {str(e)}")
                            print(f"近邻对象类型: {type(neighbor)}")
                            print(f"近邻对象内容: {neighbor}")
                            # 添加一个占位符
                            neighbors.append({
                                'element': '未知',
                                'distance': None,
                                'weight': None
                            })
                    
                    coordination_details.append({
                        'site_index': i,
                        'element': site.species_string,
                        'coordination_number': cn,
                        'neighbors': neighbors[:5]  # 只保留前5个近邻以避免数据过大
                    })
                except Exception as e:
                    print(f"计算原子 {i} 的配位环境时出错: {str(e)}")
                    print(traceback.format_exc())
            
            # 计算平均配位数
            avg_cn = sum(coordination_numbers) / len(coordination_numbers) if coordination_numbers else 0
            print(f"平均配位数: {avg_cn}")
            
            return {
                'min_distance': min_distance,  # 最短原子距离 (Å)
                'avg_coordination_number': avg_cn,  # 平均配位数
                'coordination_details': coordination_details[:5]  # 只返回前5个原子的配位细节
            }
        except Exception as e:
            print(f"计算原子距离时出错: {str(e)}")
            print(traceback.format_exc())
            return {
                'min_distance': None,
                'avg_coordination_number': None,
                'error': str(e)
            }
    
    # 新增方法: 估算表面能
    def _estimate_surface_energies(self, structure):
        """
        估算材料的表面能
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 表面能数据
        """
        try:
            print("开始估算表面能")
            # 获取最常见的晶面
            # 注意：真实的表面能计算需要DFT，这里只是生成晶面并提供几何信息
            
            # 使用SpacegroupAnalyzer获取常见的Miller指数
            spg = SpacegroupAnalyzer(structure)
            conventional_structure = spg.get_conventional_standard_structure()
            
            # 常见的低指数晶面
            miller_indices = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0), (2, 1, 1)]
            
            surface_data = []
            for miller_index in miller_indices:
                try:
                    print(f"处理晶面 {miller_index}")
                    # 生成表面切片
                    slab_gen = SlabGenerator(
                        conventional_structure, 
                        miller_index, 
                        min_slab_size=10, 
                        min_vacuum_size=15,
                        center_slab=True
                    )
                    slabs = slab_gen.get_slabs()
                    
                    if slabs:
                        slab = slabs[0]
                        # 计算表面积
                        surface_area = slab.surface_area
                        
                        # 计算表面原子数 - 修复surface_properties可能不存在的问题
                        try:
                            n_surface_atoms = len([site for site in slab if hasattr(site, 'surface_properties') and site.surface_properties])
                        except Exception as e:
                            print(f"计算表面原子数时出错: {str(e)}")
                            # 使用备用方法估算表面原子数
                            n_surface_atoms = int(len(slab) * 0.2)  # 假设约20%的原子在表面
                        
                        # 添加到结果中
                        surface_data.append({
                            'miller_index': miller_index,
                            'surface_area': surface_area,  # 表面积 (Å²)
                            'n_atoms': len(slab),  # 切片中的原子数
                            'n_surface_atoms': n_surface_atoms,  # 表面原子数
                            'thickness': slab.thickness  # 切片厚度 (Å)
                        })
                        print(f"晶面 {miller_index} 处理完成: 表面积={surface_area:.2f} Å², 原子数={len(slab)}, 表面原子数={n_surface_atoms}")
                except Exception as e:
                    print(f"生成晶面 {miller_index} 时出错: {str(e)}")
                    print(traceback.format_exc())
                    continue
            
            if not surface_data:
                print("没有成功生成任何晶面数据")
                return {
                    'error': "无法生成晶面数据",
                    'note': '注意：真实的表面能计算需要DFT，这里只提供几何信息'
                }
            
            return {
                'surfaces': surface_data,
                'note': '注意：真实的表面能计算需要DFT，这里只提供几何信息'
            }
        except Exception as e:
            print(f"估算表面能时出错: {str(e)}")
            print(traceback.format_exc())
            return {
                'error': str(e),
                'note': '无法估算表面能'
            }
    
    # 新增方法: 生成VASP输入文件
    def _generate_vasp_input(self, structure):
        """
        为结构生成VASP输入文件
        
        参数:
        structure (Structure): 材料结构
        
        返回:
        dict: 包含VASP输入文件内容的字典
        """
        try:
            # 生成POSCAR
            poscar = Poscar(structure)
            poscar_str = str(poscar)
            
            # 生成默认的INCAR
            incar_dict = {
                'SYSTEM': structure.composition.reduced_formula,
                'ENCUT': 520,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'NSW': 100,
                'ISIF': 3,
                'EDIFF': 1E-6,
                'EDIFFG': -0.01,
                'PREC': 'Accurate',
                'LREAL': 'Auto',
                'ALGO': 'Normal'
            }
            incar = Incar(incar_dict)
            incar_str = str(incar)
            
            # 生成默认的KPOINTS (自动网格)
            kpoints = Kpoints.automatic_density(structure, 1000)
            kpoints_str = str(kpoints)
            
            # 返回VASP输入文件内容
            return {
                'POSCAR': poscar_str,
                'INCAR': incar_str,
                'KPOINTS': kpoints_str,
                'note': '注意：POTCAR文件需要单独生成，因为它们受版权保护'
            }
        except Exception as e:
            print(f"生成VASP输入文件时出错: {str(e)}")
            return {
                'error': str(e),
                'note': '无法生成VASP输入文件'
            } 