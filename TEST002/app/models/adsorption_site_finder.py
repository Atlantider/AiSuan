from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io
from typing import List, Dict, Tuple, Optional, Any

class AdsorptionSiteFinder:
    """用于查找和可视化吸附位点的类"""
    
    def __init__(self):
        """初始化吸附位点查找器"""
        self.adsorbate_molecules = {
            'H': Molecule(['H'], [[0, 0, 0]]),
            'O': Molecule(['O'], [[0, 0, 0]]),
            'N': Molecule(['N'], [[0, 0, 0]]),
            'C': Molecule(['C'], [[0, 0, 0]]),
            'CO': Molecule(['C', 'O'], [[0, 0, 0], [0, 0, 1.2]]),
            'CO2': Molecule(['C', 'O', 'O'], [[0, 0, 0], [0, 0, 1.2], [0, 0, -1.2]]),
            'OH': Molecule(['O', 'H'], [[0, 0, 0], [0, 0, 0.96]]),
            'H2O': Molecule(['O', 'H', 'H'], [[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])
        }
    
    def generate_slab(self, structure: Structure, miller_index: Tuple[int, int, int], 
                     min_slab_size: float, min_vacuum_size: float, 
                     center_slab: bool = True) -> Structure:
        """
        生成表面板层
        
        Args:
            structure: 体相结构
            miller_index: 米勒指数
            min_slab_size: 最小板层厚度
            min_vacuum_size: 最小真空层厚度
            center_slab: 是否将板层居中
            
        Returns:
            生成的板层结构
        """
        # 使用pymatgen的SlabGenerator生成表面板层
        sg = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size, 
                          center_slab=center_slab, in_unit_planes=True)
        slabs = sg.get_slabs()
        
        if not slabs:
            raise ValueError(f"无法为米勒指数 {miller_index} 生成有效的表面板层")
        
        # 返回第一个板层
        return slabs[0]
    
    def find_adsorption_sites(self, slab: Structure) -> Dict[str, List[List[float]]]:
        """
        查找吸附位点
        
        Args:
            slab: 板层结构
            
        Returns:
            包含不同类型吸附位点的字典
        """
        # 使用pymatgen的AdsorbateSiteFinder查找吸附位点
        asf = AdsorbateSiteFinder(slab)
        
        # 查找所有类型的吸附位点
        sites = asf.find_adsorption_sites(distance=2.0)
        
        # 如果没有找到任何位点，创建一个默认位点
        if not sites or all(len(sites_list) == 0 for sites_list in sites.values()):
            # 找到最上层原子的z坐标
            max_z = max([site.coords[2] for site in slab])
            default_site = [slab.lattice.a/2, slab.lattice.b/2, max_z + 2.0]
            
            # 创建一个包含默认位点的字典
            sites = {'top': [default_site]}
        
        return sites
    
    def add_adsorbate(self, slab: Structure, species: str, 
                     site: List[float]) -> Structure:
        """
        在指定位点添加吸附物种
        
        Args:
            slab: 板层结构
            species: 吸附物种
            site: 吸附位点坐标
            
        Returns:
            添加吸附物种后的结构
        """
        # 获取吸附物种分子
        if species not in self.adsorbate_molecules:
            raise ValueError(f"不支持的吸附物种: {species}")
        
        molecule = self.adsorbate_molecules[species]
        
        # 使用AdsorbateSiteFinder添加吸附物种
        asf = AdsorbateSiteFinder(slab)
        ads_slab = asf.add_adsorbate(molecule, site)
        
        return ads_slab
    
    def generate_adsorption_sites_image(self, slab: Structure, 
                                      sites: Dict[str, List[List[float]]], 
                                      view_direction: str = 'c') -> bytes:
        """
        生成吸附位点图像
        
        Args:
            slab: 板层结构
            sites: 吸附位点
            view_direction: 视图方向 ('a', 'b', 或 'c')
            
        Returns:
            图像数据
        """
        # 创建图像
        fig = Figure(figsize=(8, 6), dpi=100)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        # 设置视图方向
        if view_direction == 'a':
            x_idx, y_idx = 1, 2  # b-c平面
        elif view_direction == 'b':
            x_idx, y_idx = 0, 2  # a-c平面
        else:  # 默认为c方向
            x_idx, y_idx = 0, 1  # a-b平面
        
        # 绘制板层原子
        elements = set([site.specie.symbol for site in slab])
        element_colors = {}
        
        # 为每个元素分配颜色
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray']
        for i, element in enumerate(elements):
            element_colors[element] = colors[i % len(colors)]
        
        # 绘制板层原子
        for site in slab:
            element = site.specie.symbol
            coords = site.coords
            ax.scatter(coords[x_idx], coords[y_idx], c=element_colors[element], 
                      s=100, label=element, alpha=0.7)
        
        # 绘制吸附位点
        site_colors = {'top': 'red', 'bridge': 'green', 'hollow': 'blue', 'ontop': 'red'}
        site_markers = {'top': 'o', 'bridge': 's', 'hollow': '^', 'ontop': 'o'}
        
        for site_type, site_list in sites.items():
            color = site_colors.get(site_type, 'black')
            marker = site_markers.get(site_type, 'x')
            
            for site in site_list:
                ax.scatter(site[x_idx], site[y_idx], c=color, marker=marker, 
                          s=80, label=f"{site_type} site", alpha=0.5)
        
        # 移除重复的图例
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='best')
        
        # 设置标题和轴标签
        ax.set_title(f'吸附位点 ({view_direction}方向视图)')
        ax.set_xlabel(f'{["b", "a", "a"][["a", "b", "c"].index(view_direction)]} (Å)')
        ax.set_ylabel(f'{["c", "c", "b"][["a", "b", "c"].index(view_direction)]} (Å)')
        
        # 调整布局
        fig.tight_layout()
        
        # 将图像保存到内存中
        buf = io.BytesIO()
        canvas.print_png(buf)
        buf.seek(0)
        
        return buf.getvalue()
    
    def generate_adsorption_structure_image(self, structure: Structure, 
                                          view_direction: str = 'c') -> bytes:
        """
        生成吸附结构图像
        
        Args:
            structure: 吸附结构
            view_direction: 视图方向 ('a', 'b', 或 'c')
            
        Returns:
            图像数据
        """
        # 创建图像
        fig = Figure(figsize=(8, 6), dpi=100)
        canvas = FigureCanvas(fig)
        ax = fig.add_subplot(111)
        
        # 设置视图方向
        if view_direction == 'a':
            x_idx, y_idx = 1, 2  # b-c平面
        elif view_direction == 'b':
            x_idx, y_idx = 0, 2  # a-c平面
        else:  # 默认为c方向
            x_idx, y_idx = 0, 1  # a-b平面
        
        # 绘制结构原子
        elements = set([site.specie.symbol for site in structure])
        element_colors = {}
        
        # 为每个元素分配颜色
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray']
        for i, element in enumerate(elements):
            element_colors[element] = colors[i % len(colors)]
        
        # 找到最大的z坐标，用于区分吸附物种和表面
        max_slab_z = np.percentile([site.coords[2] for site in structure], 90)
        
        # 绘制结构原子
        for site in structure:
            element = site.specie.symbol
            coords = site.coords
            
            # 区分吸附物种和表面原子
            if coords[2] > max_slab_z:
                # 吸附物种
                ax.scatter(coords[x_idx], coords[y_idx], c=element_colors[element], 
                          s=120, edgecolors='black', label=f"{element} (吸附物种)", alpha=0.9)
            else:
                # 表面原子
                ax.scatter(coords[x_idx], coords[y_idx], c=element_colors[element], 
                          s=100, label=element, alpha=0.7)
        
        # 移除重复的图例
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='best')
        
        # 设置标题和轴标签
        ax.set_title(f'吸附结构 ({view_direction}方向视图)')
        ax.set_xlabel(f'{["b", "a", "a"][["a", "b", "c"].index(view_direction)]} (Å)')
        ax.set_ylabel(f'{["c", "c", "b"][["a", "b", "c"].index(view_direction)]} (Å)')
        
        # 调整布局
        fig.tight_layout()
        
        # 将图像保存到内存中
        buf = io.BytesIO()
        canvas.print_png(buf)
        buf.seek(0)
        
        return buf.getvalue() 