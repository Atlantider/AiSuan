"""
分子计算器模块

根据物理化学参数（体积、浓度、比例等）自动计算模拟所需的分子数量，
无需用户在Excel中手动计算。
"""

import logging
import math
from typing import Dict, List, Tuple, Optional, Union, Any
import numpy as np

class MoleculeCalculator:
    """
    分子数量计算器
    
    根据物理化学参数自动计算模拟中需要的分子数量
    """
    
    # 阿伏伽德罗常数
    AVOGADRO = 6.02214076e23  # 单位：mol^-1
    
    def __init__(self, logger=None):
        """
        初始化分子计算器
        
        Args:
            logger: 日志记录器实例
        """
        self.logger = logger or logging.getLogger(__name__)
    
    def calculate_molecule_counts(self, 
                               volume: float, 
                               concentration: float,
                               solvent_ratio: float = 0.0,
                               cation_anion_ratio: float = 1.0,
                               volume_unit: str = 'nm3', 
                               concentration_unit: str = 'mol/L') -> Dict[str, int]:
        """
        计算所需的离子和溶剂数量
        
        Args:
            volume: 模拟盒子的体积
            concentration: 离子的摩尔浓度
            solvent_ratio: 溶剂与离子的比例 (溶剂分子数/离子对数)
            cation_anion_ratio: 阳离子与阴离子的比例 (默认为1:1)
            volume_unit: 体积单位，可选 'nm3' (立方纳米) 或 'angstrom3' (立方埃)
            concentration_unit: 浓度单位，可选 'mol/L'
            
        Returns:
            包含阳离子、阴离子和溶剂数量的字典
        """
        self.logger.info(f"计算分子数量：体积={volume}{volume_unit}，浓度={concentration}{concentration_unit}")
        
        # 转换体积为立方米
        volume_m3 = self._convert_volume_to_m3(volume, volume_unit)
        
        # 转换浓度为mol/m3
        concentration_mol_m3 = self._convert_concentration(concentration, concentration_unit)
        
        # 计算总离子对数量
        ion_pairs = concentration_mol_m3 * volume_m3 * self.AVOGADRO
        
        # 向上取整确保有足够的离子对
        ion_pairs_int = math.ceil(ion_pairs)
        
        # 考虑阳离子与阴离子的比例
        if cation_anion_ratio >= 1.0:
            cation_count = ion_pairs_int
            anion_count = math.ceil(ion_pairs_int / cation_anion_ratio)
        else:
            anion_count = ion_pairs_int
            cation_count = math.ceil(ion_pairs_int * cation_anion_ratio)
        
        # 计算溶剂分子数量
        solvent_count = math.ceil(ion_pairs_int * solvent_ratio)
        
        result = {
            'cation': cation_count,
            'anion': anion_count,
            'solvent': solvent_count,
            'total': cation_count + anion_count + solvent_count
        }
        
        self.logger.info(f"计算结果：阳离子={result['cation']}，阴离子={result['anion']}，溶剂={result['solvent']}")
        return result
    
    def calculate_from_system_definition(self, 
                                      box_size: Tuple[float, float, float], 
                                      cation_concentration: float, 
                                      solvent_data: Optional[Dict[str, float]] = None, 
                                      cation_types: int = 1,
                                      anion_types: int = 1,
                                      length_unit: str = 'angstrom',
                                      concentration_unit: str = 'mol/L') -> Dict[str, Dict[str, int]]:
        """
        从模拟系统定义计算所有分子数量
        
        Args:
            box_size: 模拟盒子的尺寸 (x, y, z)
            cation_concentration: 阳离子的摩尔浓度
            solvent_data: 溶剂数据，格式为 {溶剂名称: 溶剂比例}
            cation_types: 阳离子类型数量
            anion_types: 阴离子类型数量
            length_unit: 长度单位，可选 'angstrom' 或 'nm'
            concentration_unit: 浓度单位，可选 'mol/L'
            
        Returns:
            包含所有分子类型数量的嵌套字典
        """
        self.logger.info(f"从系统定义计算分子数量: 盒子尺寸={box_size}{length_unit}")
        
        # 计算盒子体积
        volume = box_size[0] * box_size[1] * box_size[2]
        volume_unit = 'angstrom3' if length_unit == 'angstrom' else 'nm3'
        
        # 计算总离子和溶剂数量
        total_solvent_ratio = sum(solvent_data.values()) if solvent_data else 0.0
        molecule_counts = self.calculate_molecule_counts(
            volume=volume, 
            concentration=cation_concentration,
            solvent_ratio=total_solvent_ratio,
            volume_unit=volume_unit,
            concentration_unit=concentration_unit
        )
        
        # 计算每种离子的详细数量
        result = {
            'cations': {},
            'anions': {},
            'solvents': {}
        }
        
        # 分配阳离子
        if cation_types == 1:
            result['cations']['cation1'] = molecule_counts['cation']
        else:
            # 平均分配阳离子，或者可以实现更复杂的分配逻辑
            base_count = molecule_counts['cation'] // cation_types
            remainder = molecule_counts['cation'] % cation_types
            
            for i in range(1, cation_types + 1):
                key = f'cation{i}'
                if i <= remainder:
                    result['cations'][key] = base_count + 1
                else:
                    result['cations'][key] = base_count
        
        # 分配阴离子
        if anion_types == 1:
            result['anions']['anion1'] = molecule_counts['anion']
        else:
            # 平均分配阴离子
            base_count = molecule_counts['anion'] // anion_types
            remainder = molecule_counts['anion'] % anion_types
            
            for i in range(1, anion_types + 1):
                key = f'anion{i}'
                if i <= remainder:
                    result['anions'][key] = base_count + 1
                else:
                    result['anions'][key] = base_count
        
        # 分配溶剂
        if solvent_data:
            total_ratio = sum(solvent_data.values())
            
            for solvent_name, ratio in solvent_data.items():
                if total_ratio > 0:
                    solvent_count = math.ceil(molecule_counts['solvent'] * (ratio / total_ratio))
                else:
                    solvent_count = 0
                
                result['solvents'][solvent_name] = solvent_count
        
        # 记录总数
        result['total'] = sum(result['cations'].values()) + sum(result['anions'].values()) + sum(result['solvents'].values())
        
        self.logger.info(f"计算完成，总分子数: {result['total']}")
        return result
    
    def calculate_custom_composition(self, 
                                  box_size: Tuple[float, float, float],
                                  components: Dict[str, Dict[str, Any]],
                                  length_unit: str = 'angstrom') -> Dict[str, Dict[str, int]]:
        """
        计算自定义组分的分子数量
        
        Args:
            box_size: 模拟盒子的尺寸 (x, y, z)
            components: 组分数据，格式为 {组分名称: {
                'type': '类型（cation/anion/solvent）', 
                'concentration': 浓度, 
                'concentration_unit': 浓度单位,
                'ratio': 比例（可选）}
            }
            length_unit: 长度单位，可选 'angstrom' 或 'nm'
            
        Returns:
            包含所有分子类型数量的嵌套字典
        """
        self.logger.info(f"计算自定义组分的分子数量: 盒子尺寸={box_size}{length_unit}")
        
        # 计算盒子体积
        volume = box_size[0] * box_size[1] * box_size[2]
        volume_unit = 'angstrom3' if length_unit == 'angstrom' else 'nm3'
        
        result = {
            'cations': {},
            'anions': {},
            'solvents': {}
        }
        
        # 首先处理基于浓度的组分
        for name, data in components.items():
            if 'concentration' in data:
                concentration = data['concentration']
                concentration_unit = data.get('concentration_unit', 'mol/L')
                
                # 转换体积为立方米
                volume_m3 = self._convert_volume_to_m3(volume, volume_unit)
                
                # 转换浓度为mol/m3
                concentration_mol_m3 = self._convert_concentration(concentration, concentration_unit)
                
                # 计算分子数量
                molecule_count = math.ceil(concentration_mol_m3 * volume_m3 * self.AVOGADRO)
                
                # 根据类型添加到结果中
                comp_type = data.get('type', 'solvent').lower()
                if comp_type == 'cation':
                    result['cations'][name] = molecule_count
                elif comp_type == 'anion':
                    result['anions'][name] = molecule_count
                else:
                    result['solvents'][name] = molecule_count
        
        # 然后处理基于比例的组分
        ratio_components = {name: data for name, data in components.items() if 'ratio' in data}
        
        if ratio_components:
            # 计算总离子数
            total_ions = sum(result['cations'].values()) + sum(result['anions'].values())
            
            for name, data in ratio_components.items():
                ratio = data['ratio']
                comp_type = data.get('type', 'solvent').lower()
                
                if comp_type == 'solvent':
                    # 溶剂比例通常是相对于离子对数量
                    molecule_count = math.ceil(total_ions * ratio / 2)  # 除以2是因为一般以离子对为单位
                    result['solvents'][name] = molecule_count
                else:
                    # 对于阳离子或阴离子，比例通常是相对于总离子数
                    molecule_count = math.ceil(total_ions * ratio)
                    if comp_type == 'cation':
                        result['cations'][name] = molecule_count
                    else:  # anion
                        result['anions'][name] = molecule_count
        
        # 计算总数
        result['total'] = sum(result['cations'].values()) + sum(result['anions'].values()) + sum(result['solvents'].values())
        
        self.logger.info(f"自定义组分计算完成，总分子数: {result['total']}")
        return result
    
    def calculate_from_primary_cation(self,
                                   box_size: float,
                                   primary_cation_concentration: float,
                                   component_ratios: Dict[str, float],
                                   length_unit: str = 'angstrom',
                                   concentration_unit: str = 'mol/L') -> Dict[str, Dict[str, int]]:
        """
        基于主阳离子和比例计算所有分子数量
        
        以第一个阳离子为基准(ratio=1)，其他所有组分按照相对比例计算
        
        Args:
            box_size: 立方体盒子的边长
            primary_cation_concentration: 主阳离子的浓度
            component_ratios: 各组分相对于主阳离子的比例，格式为 {组分名称: 比例值}
                              组分名称应该包含类型前缀：cation、anion或sol
            length_unit: 长度单位，可选 'angstrom' 或 'nm'
            concentration_unit: 浓度单位，可选 'mol/L'
            
        Returns:
            包含所有分子类型数量的嵌套字典
        """
        self.logger.info(f"基于主阳离子计算分子数量: 盒子尺寸={box_size}{length_unit}, 浓度={primary_cation_concentration}{concentration_unit}")
        
        # 计算立方体盒子体积
        volume = box_size ** 3
        volume_unit = 'angstrom3' if length_unit == 'angstrom' else 'nm3'
        
        # 转换体积为立方米
        volume_m3 = self._convert_volume_to_m3(volume, volume_unit)
        
        # 转换浓度为mol/m3
        concentration_mol_m3 = self._convert_concentration(primary_cation_concentration, concentration_unit)
        
        # 计算主阳离子数量
        primary_cation_count = math.ceil(concentration_mol_m3 * volume_m3 * self.AVOGADRO)
        self.logger.info(f"主阳离子数量: {primary_cation_count}")
        
        # 分类存储计算结果
        result = {
            'cations': {},
            'anions': {},
            'solvents': {}
        }
        
        # 添加主阳离子
        result['cations']['cation1'] = primary_cation_count
        
        # 根据比例计算其他组分数量
        for component_name, ratio in component_ratios.items():
            # 根据前缀确定组分类型
            if component_name.startswith('cation'):
                if component_name == 'cation1':  # 主阳离子已经处理
                    continue
                # 计算数量并添加到阳离子列表
                count = math.ceil(primary_cation_count * ratio)
                result['cations'][component_name] = count
                self.logger.info(f"{component_name} 数量: {count} (比例: {ratio})")
            elif component_name.startswith('anion'):
                # 计算数量并添加到阴离子列表
                count = math.ceil(primary_cation_count * ratio)
                result['anions'][component_name] = count
                self.logger.info(f"{component_name} 数量: {count} (比例: {ratio})")
            elif component_name.startswith('sol'):
                # 计算数量并添加到溶剂列表
                count = math.ceil(primary_cation_count * ratio)
                result['solvents'][component_name] = count
                self.logger.info(f"{component_name} 数量: {count} (比例: {ratio})")
            else:
                self.logger.warning(f"未知组分类型: {component_name}")
        
        # 计算总数
        result['total'] = sum(result['cations'].values()) + sum(result['anions'].values()) + sum(result['solvents'].values())
        
        self.logger.info(f"基于主阳离子的计算完成，总分子数: {result['total']}")
        return result
    
    def estimate_density(self, 
                       molecule_counts: Dict[str, Dict[str, int]], 
                       molecular_weights: Dict[str, float],
                       box_size: Tuple[float, float, float],
                       length_unit: str = 'angstrom') -> float:
        """
        估计模拟体系的密度
        
        Args:
            molecule_counts: 分子数量字典
            molecular_weights: 各分子的分子量 (g/mol)
            box_size: 模拟盒子的尺寸 (x, y, z)
            length_unit: 长度单位，可选 'angstrom' 或 'nm'
            
        Returns:
            估计的密度 (g/cm³)
        """
        # 计算盒子体积
        volume = box_size[0] * box_size[1] * box_size[2]
        
        # 转换体积为立方厘米
        if length_unit == 'angstrom':
            volume_cm3 = volume * 1e-24  # 立方埃 -> 立方厘米
        else:  # 'nm'
            volume_cm3 = volume * 1e-21  # 立方纳米 -> 立方厘米
        
        # 计算总质量 (g)
        total_mass = 0.0
        
        # 处理阳离子
        for name, count in molecule_counts.get('cations', {}).items():
            if name in molecular_weights:
                total_mass += count * molecular_weights[name] / self.AVOGADRO
        
        # 处理阴离子
        for name, count in molecule_counts.get('anions', {}).items():
            if name in molecular_weights:
                total_mass += count * molecular_weights[name] / self.AVOGADRO
        
        # 处理溶剂
        for name, count in molecule_counts.get('solvents', {}).items():
            if name in molecular_weights:
                total_mass += count * molecular_weights[name] / self.AVOGADRO
        
        # 计算密度 (g/cm³)
        if volume_cm3 > 0:
            density = total_mass / volume_cm3
        else:
            density = 0.0
        
        self.logger.info(f"估计系统密度: {density:.4f} g/cm³")
        return density
    
    def _convert_volume_to_m3(self, volume: float, unit: str) -> float:
        """
        将体积转换为立方米
        
        Args:
            volume: 体积值
            unit: 体积单位
            
        Returns:
            立方米单位的体积
        """
        if unit == 'nm3':
            return volume * 1e-27  # 立方纳米 -> 立方米
        elif unit == 'angstrom3':
            return volume * 1e-30  # 立方埃 -> 立方米
        elif unit == 'm3':
            return volume
        else:
            self.logger.warning(f"未知的体积单位: {unit}，假设为立方米")
            return volume
    
    def _convert_concentration(self, concentration: float, unit: str) -> float:
        """
        将浓度转换为mol/m3
        
        Args:
            concentration: 浓度值
            unit: 浓度单位
            
        Returns:
            mol/m3单位的浓度
        """
        if unit == 'mol/L':
            return concentration * 1000  # mol/L -> mol/m3
        elif unit == 'mol/m3':
            return concentration
        else:
            self.logger.warning(f"未知的浓度单位: {unit}，假设为mol/m3")
            return concentration 