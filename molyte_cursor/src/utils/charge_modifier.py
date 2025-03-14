"""
电荷修改器模块，用于修改LAMMPS文件中的电荷信息
"""
from pathlib import Path
from .logger import Logger

class ChargeModifier:
    """电荷修改器类"""
    
    def __init__(self):
        """初始化电荷修改器"""
        self.logger = Logger().get_logger()
    
    def read_file(self, file_path):
        """读取文件内容
        
        Args:
            file_path: 文件路径
            
        Returns:
            文件内容行列表
        """
        with open(file_path, 'r') as file:
            return file.readlines()
    
    def write_file(self, file_path, contents):
        """写入文件内容
        
        Args:
            file_path: 文件路径
            contents: 要写入的内容行列表
        """
        with open(file_path, 'w') as file:
            file.writelines(contents)
    
    def extract_charges(self, chg_contents):
        """从电荷文件中提取电荷值
        
        Args:
            chg_contents: 电荷文件内容行列表
            
        Returns:
            电荷值列表
        """
        charges = []
        for line in chg_contents:
            parts = line.split()
            if parts:
                try:
                    charge = float(parts[-1])
                    charges.append(charge)
                except ValueError:
                    continue
        return charges
    
    def find_atoms_section(self, lmp_contents):
        """在LAMMPS文件中查找原子部分
        
        Args:
            lmp_contents: LAMMPS文件内容行列表
            
        Returns:
            (原子部分开始索引, 原子部分结束索引)
        """
        atom_section_start = None
        atom_section_end = None
        
        for i, line in enumerate(lmp_contents):
            if line == "Atoms\n":
                atom_section_start = i + 2
                break
                
        if atom_section_start is not None:
            for i in range(atom_section_start, len(lmp_contents)):
                if lmp_contents[i].strip() == '':
                    atom_section_end = i
                    break
        
        return atom_section_start, atom_section_end
    
    def modify_charges_in_lmp(self, lmp_contents, atom_section_start, atom_section_end, chg_values):
        """在LAMMPS文件中修改电荷值
        
        Args:
            lmp_contents: LAMMPS文件内容行列表
            atom_section_start: 原子部分开始索引
            atom_section_end: 原子部分结束索引
            chg_values: 新的电荷值列表
            
        Returns:
            修改后的LAMMPS文件内容行列表
        """
        modified_section = []
        for line, new_charge in zip(lmp_contents[atom_section_start:atom_section_end], chg_values):
            parts = line.split()
            if parts:
                parts[3] = str(new_charge)
                modified_line = ' '.join(parts) + '\n'
                modified_section.append(modified_line)
            else:
                modified_section.append(line)
        
        lmp_contents[atom_section_start:atom_section_end] = modified_section
        return lmp_contents
    
    def modify_lmp_charges(self, lmp_file_path, chg_file_path, output_lmp_file_path):
        """修改LAMMPS文件中的电荷
        
        Args:
            lmp_file_path: LAMMPS文件路径
            chg_file_path: 电荷文件路径
            output_lmp_file_path: 输出LAMMPS文件路径
            
        Raises:
            ValueError: 当原子数量和电荷数量不匹配时
        """
        # 读取文件内容
        lmp_contents = self.read_file(lmp_file_path)
        chg_contents = self.read_file(chg_file_path)
        
        # 提取电荷值
        chg_values = self.extract_charges(chg_contents)
        
        # 查找原子部分
        atom_section_start, atom_section_end = self.find_atoms_section(lmp_contents)
        
        # 确保原子数量与电荷数量匹配
        if len(lmp_contents[atom_section_start:atom_section_end]) != len(chg_values):
            raise ValueError("Number of atoms and charges do not match.")
        
        # 修改LAMMPS文件中的电荷
        modified_lmp_contents = self.modify_charges_in_lmp(
            lmp_contents, atom_section_start, atom_section_end, chg_values
        )
        
        # 写入修改后的LAMMPS文件
        self.write_file(output_lmp_file_path, modified_lmp_contents)
        self.logger.info(f"Modified charges in {output_lmp_file_path}") 