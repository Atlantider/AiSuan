from pymatgen.core.structure import Structure, IStructure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.vis.structure_vtk import StructureVis
import numpy as np
from typing import Tuple, List, Optional
import os
import io
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import tempfile

class StructureHandler:
    def __init__(self):
        self.current_structure = None

    def load_structure(self, file_path: str) -> bool:
        """从文件加载晶体结构"""
        try:
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext == '.cif':
                parser = CifParser(file_path)
                self.current_structure = parser.get_structures()[0]
            elif file_ext in ['.vasp', '.poscar']:
                self.current_structure = Poscar.from_file(file_path).structure
            return True
        except Exception as e:
            print(f"加载结构失败: {str(e)}")
            return False

    def create_structure_from_parameters(
        self,
        lattice_params: Tuple[float, float, float],
        angles: Tuple[float, float, float],
        species: List[str],
        coords: List[List[float]]
    ) -> bool:
        """从晶格参数创建结构"""
        try:
            a, b, c = lattice_params
            alpha, beta, gamma = angles
            
            # 创建晶格矩阵
            lattice = Structure.from_parameters(a, b, c, alpha, beta, gamma)
            
            # 创建结构
            self.current_structure = Structure(
                lattice=lattice.lattice,
                species=species,
                coords=coords,
                coords_are_cartesian=False
            )
            return True
        except Exception as e:
            print(f"创建结构失败: {str(e)}")
            return False

    def create_supercell(self, scaling_matrix: Tuple[int, int, int]) -> bool:
        """创建超胞"""
        try:
            if self.current_structure is None:
                return False
            self.current_structure.make_supercell(scaling_matrix)
            return True
        except Exception as e:
            print(f"创建超胞失败: {str(e)}")
            return False

    def add_atom(self, species: str, coords: List[float], coords_are_cartesian: bool = False) -> bool:
        """添加原子"""
        try:
            if self.current_structure is None:
                return False
            self.current_structure.append(species, coords, coords_are_cartesian)
            return True
        except Exception as e:
            print(f"添加原子失败: {str(e)}")
            return False

    def delete_atom(self, index: int) -> bool:
        """删除原子"""
        try:
            if self.current_structure is None or index >= len(self.current_structure):
                return False
            self.current_structure.remove_sites([index])
            return True
        except Exception as e:
            print(f"删除原子失败: {str(e)}")
            return False

    def modify_atom_position(self, index: int, new_coords: List[float], coords_are_cartesian: bool = False) -> bool:
        """修改原子位置"""
        try:
            if self.current_structure is None or index >= len(self.current_structure):
                return False
            if coords_are_cartesian:
                self.current_structure.cart_coords[index] = new_coords
            else:
                self.current_structure.frac_coords[index] = new_coords
            return True
        except Exception as e:
            print(f"修改原子位置失败: {str(e)}")
            return False

    def create_immutable_structure(self) -> Optional[IStructure]:
        """创建不可变结构对象"""
        try:
            if self.current_structure is None:
                return None
            return IStructure.from_sites(self.current_structure.sites)
        except Exception as e:
            print(f"创建不可变结构失败: {str(e)}")
            return None

    def get_structure_dict(self) -> dict:
        """获取结构的字典表示"""
        if self.current_structure is None:
            return {}
        
        return {
            'lattice': self.current_structure.lattice.matrix.tolist(),
            'species': [str(site.specie) for site in self.current_structure],
            'coords': self.current_structure.frac_coords.tolist(),
            'cart_coords': self.current_structure.cart_coords.tolist(),
            'formula': self.current_structure.formula,
            'num_sites': len(self.current_structure)
        }
    
    def generate_structure_image(self, structure_dict: dict = None, view_direction: str = 'a') -> Optional[bytes]:
        """生成结构图像"""
        try:
            # 使用当前结构或传入的结构字典
            structure = self.current_structure
            if structure is None and structure_dict:
                # 从字典创建结构
                lattice = structure_dict['lattice']
                species = structure_dict['species']
                coords = structure_dict['coords']
                
                structure = Structure(
                    lattice=lattice,
                    species=species,
                    coords=coords,
                    coords_are_cartesian=False
                )
            
            if structure is None:
                return None
            
            # 创建图像
            fig = Figure(figsize=(8, 8), dpi=100)
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            
            # 设置视图方向
            if view_direction == 'a':
                proj_matrix = np.array([[0, 1, 0], [0, 0, 1]])  # bc平面 (垂直于a轴)
            elif view_direction == 'b':
                proj_matrix = np.array([[1, 0, 0], [0, 0, 1]])  # ac平面 (垂直于b轴)
            else:  # 'c'
                proj_matrix = np.array([[1, 0, 0], [0, 1, 0]])  # ab平面 (垂直于c轴)
            
            # 获取晶格向量
            lattice = structure.lattice.matrix
            
            # 投影晶格向量
            if view_direction == 'a':
                proj_vectors = [np.dot(proj_matrix, lattice[1]), np.dot(proj_matrix, lattice[2])]
            elif view_direction == 'b':
                proj_vectors = [np.dot(proj_matrix, lattice[0]), np.dot(proj_matrix, lattice[2])]
            else:  # 'c'
                proj_vectors = [np.dot(proj_matrix, lattice[0]), np.dot(proj_matrix, lattice[1])]
            
            # 绘制晶格
            origin = np.array([0, 0])
            ax.arrow(origin[0], origin[1], proj_vectors[0][0], proj_vectors[0][1], 
                    head_width=0.1, head_length=0.2, fc='black', ec='black', length_includes_head=True)
            ax.arrow(origin[0], origin[1], proj_vectors[1][0], proj_vectors[1][1], 
                    head_width=0.1, head_length=0.2, fc='black', ec='black', length_includes_head=True)
            
            # 绘制原子
            colors = {'H': 'white', 'C': 'black', 'N': 'blue', 'O': 'red', 'F': 'green', 
                    'Na': 'purple', 'Mg': 'pink', 'Al': 'gray', 'Si': 'orange', 'P': 'brown', 
                    'S': 'yellow', 'Cl': 'green', 'K': 'purple', 'Ca': 'gray', 'Ti': 'gray', 
                    'Cr': 'gray', 'Mn': 'gray', 'Fe': 'brown', 'Co': 'blue', 'Ni': 'gray', 
                    'Cu': 'orange', 'Zn': 'gray', 'Ga': 'gray', 'Ge': 'gray', 'As': 'gray', 
                    'Se': 'gray', 'Br': 'brown', 'Rb': 'gray', 'Sr': 'gray', 'Y': 'gray', 
                    'Zr': 'gray', 'Nb': 'gray', 'Mo': 'gray', 'Tc': 'gray', 'Ru': 'gray', 
                    'Rh': 'gray', 'Pd': 'gray', 'Ag': 'gray', 'Cd': 'gray', 'In': 'gray', 
                    'Sn': 'gray', 'Sb': 'gray', 'Te': 'gray', 'I': 'purple', 'Cs': 'gray', 
                    'Ba': 'gray', 'La': 'gray', 'Ce': 'gray', 'Pr': 'gray', 'Nd': 'gray', 
                    'Pm': 'gray', 'Sm': 'gray', 'Eu': 'gray', 'Gd': 'gray', 'Tb': 'gray', 
                    'Dy': 'gray', 'Ho': 'gray', 'Er': 'gray', 'Tm': 'gray', 'Yb': 'gray', 
                    'Lu': 'gray', 'Hf': 'gray', 'Ta': 'gray', 'W': 'gray', 'Re': 'gray', 
                    'Os': 'gray', 'Ir': 'gray', 'Pt': 'gray', 'Au': 'yellow', 'Hg': 'gray', 
                    'Tl': 'gray', 'Pb': 'gray', 'Bi': 'gray', 'Po': 'gray', 'At': 'gray', 
                    'Rn': 'gray', 'Fr': 'gray', 'Ra': 'gray', 'Ac': 'gray', 'Th': 'gray', 
                    'Pa': 'gray', 'U': 'gray'}
            
            # 默认颜色
            default_color = 'gray'
            
            # 获取分数坐标
            frac_coords = structure.frac_coords
            
            # 投影原子坐标
            for i, site in enumerate(structure):
                # 获取元素符号
                element = str(site.specie.symbol)
                
                # 获取颜色
                color = colors.get(element, default_color)
                
                # 获取分数坐标
                fcoord = frac_coords[i]
                
                # 计算投影坐标
                if view_direction == 'a':
                    proj_coord = fcoord[1] * proj_vectors[0] + fcoord[2] * proj_vectors[1]
                elif view_direction == 'b':
                    proj_coord = fcoord[0] * proj_vectors[0] + fcoord[2] * proj_vectors[1]
                else:  # 'c'
                    proj_coord = fcoord[0] * proj_vectors[0] + fcoord[1] * proj_vectors[1]
                
                # 绘制原子
                ax.scatter(proj_coord[0], proj_coord[1], color=color, s=100, edgecolors='black', zorder=2)
                
                # 添加元素标签
                ax.text(proj_coord[0], proj_coord[1], element, fontsize=8, 
                       ha='center', va='center', color='white', zorder=3)
            
            # 设置图像属性
            ax.set_aspect('equal')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_title(f'结构视图 ({view_direction}轴方向)')
            
            # 调整坐标轴范围
            ax.autoscale_view()
            
            # 保存图像到内存
            buf = io.BytesIO()
            canvas.print_png(buf)
            buf.seek(0)
            
            return buf.getvalue()
            
        except Exception as e:
            print(f"生成结构图像失败: {str(e)}")
            return None 

    def process_uploaded_file(self, file):
        """处理上传的结构文件"""
        try:
            # 保存文件
            filename = file.filename
            file_path = os.path.join('uploads', filename)
            file.save(file_path)
            
            # 加载结构
            self.load_structure(file_path)
            
            # 返回结构数据和文件路径
            structure_data = self.get_structure_dict()
            structure_data['file_path'] = file_path
            
            return structure_data
        except Exception as e:
            raise ValueError(f"处理上传文件时出错: {str(e)}") 