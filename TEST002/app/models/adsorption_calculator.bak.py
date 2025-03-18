import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.sites import PeriodicSite
import os
import json

class AdsorptionCalculator:
    """材料表面吸附能计算器"""
    
    def __init__(self):
        """初始化计算器"""
        # 预设的吸附物种能量（eV），实际应用中应从DFT计算获取
        self.species_energies = {
            'H': -3.39,
            'O': -4.94,
            'N': -8.32,
            'C': -7.28,
            'CO': -14.80,
            'CO2': -22.96,
            'OH': -9.87,
            'H2O': -14.22
        }
    
    def calculate(self, data):
        """
        计算吸附能
        
        参数:
        data (dict): 包含计算参数的字典
            - structure_file_path: 结构文件路径
            - species: 吸附物种
            - miller_index: Miller指数，如"1,1,1"
            - slab_thickness: 板层厚度（原子层数）
            - vacuum_thickness: 真空层厚度（埃）
            - adsorption_site: 吸附位点类型（top/bridge/hollow）
            
        返回:
        dict: 计算结果
        """
        try:
            # 解析参数
            structure_path = data.get('structure_file_path')
            species = data.get('species')
            miller_index = [int(i) for i in data.get('miller_index', '1,1,1').split(',')]
            slab_thickness = int(data.get('slab_thickness', 4))
            vacuum_thickness = float(data.get('vacuum_thickness', 15.0))
            adsorption_site_type = data.get('adsorption_site', 'top')
            
            # 加载结构
            structure = Structure.from_file(structure_path)
            
            # 生成表面板层
            slabs = self._generate_slabs(structure, miller_index, slab_thickness, vacuum_thickness)
            
            if not slabs:
                return {"error": "无法生成有效的表面板层"}
            
            # 选择第一个板层（实际应用中可能需要更复杂的选择逻辑）
            slab = slabs[0]
            
            # 找到吸附位点
            sites = self._find_adsorption_sites(slab, adsorption_site_type)
            
            if not sites:
                return {"error": f"在所选表面上找不到{adsorption_site_type}类型的吸附位点"}
            
            # 选择第一个吸附位点（实际应用中可能需要更复杂的选择逻辑）
            site = sites[0]
            
            # 计算吸附能
            adsorption_energy = self._calculate_adsorption_energy(slab, site, species)
            
            # 准备结果
            result = {
                "adsorption_energy": adsorption_energy,
                "adsorption_energy_formatted": f"{adsorption_energy:.3f} eV",
                "surface": f"{miller_index[0]}{miller_index[1]}{miller_index[2]}",
                "species": species,
                "site_type": adsorption_site_type,
                "site_position": site.coords.tolist()
            }
            
            return result
            
        except Exception as e:
            return {"error": f"计算过程中出错: {str(e)}"}
    
    def _generate_slabs(self, structure, miller_index, thickness, vacuum):
        """生成表面板层"""
        # 使用pymatgen的SlabGenerator生成表面板层
        sg = SlabGenerator(structure, miller_index, thickness, vacuum, center_slab=True)
        slabs = sg.get_slabs()
        return slabs
    
    def _find_adsorption_sites(self, slab, site_type):
        """找到吸附位点"""
        # 使用pymatgen的AdsorbateSiteFinder找到吸附位点
        asf = AdsorbateSiteFinder(slab)
        
        # 根据位点类型获取位点
        if site_type == 'top':
            sites = asf.find_adsorption_sites(distance=2.0)['top']
        elif site_type == 'bridge':
            sites = asf.find_adsorption_sites(distance=2.0)['bridge']
        elif site_type == 'hollow':
            sites = asf.find_adsorption_sites(distance=2.0)['hollow']
        else:
            sites = []
        
        return sites
    
    def _calculate_adsorption_energy(self, slab, site, species):
        """
        计算吸附能
        
        这里使用简化模型：E_ads = E_slab+adsorbate - E_slab - E_adsorbate
        实际应用中应使用DFT计算获取更准确的能量
        """
        # 在实际应用中，这些能量应该通过DFT计算获得
        # 这里使用模拟值进行演示
        
        # 获取物种能量
        species_energy = self.species_energies.get(species, 0)
        
        # 模拟板层能量（实际应通过DFT计算）
        # 这里简单地基于原子数量估算
        slab_energy = -5.0 * len(slab)
        
        # 模拟吸附后的能量
        # 这里使用一个简单的模型，考虑吸附位点与表面的距离
        # 实际应用中应通过DFT计算获取
        distance_factor = np.min([np.linalg.norm(site.coords - atom.coords) for atom in slab])
        interaction_energy = -2.0 / (1.0 + 0.1 * distance_factor)
        
        slab_adsorbate_energy = slab_energy + species_energy + interaction_energy
        
        # 计算吸附能
        adsorption_energy = slab_adsorbate_energy - slab_energy - species_energy
        
        return adsorption_energy 