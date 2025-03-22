"""
电荷修改器模块，用于修改LAMMPS文件中的电荷信息
"""
from pathlib import Path
from .logger import Logger
import os
import logging

logger = logging.getLogger(__name__)

# 添加模块级别的函数
def modify_lmp_charges(lmp_file_path, chg_file_path, output_lmp_file_path):
    """
    修改LAMMPS文件中的电荷
    
    Args:
        lmp_file_path: LAMMPS文件路径
        chg_file_path: 电荷文件路径
        output_lmp_file_path: 输出LAMMPS文件路径
        
    Raises:
        FileNotFoundError: 当文件不存在时
        ValueError: 当原子数量与电荷值数量不匹配时
    """
    modifier = ChargeModifier()
    return modifier.modify_lmp_charges(lmp_file_path, chg_file_path, output_lmp_file_path)

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
            if line.strip() == "Atoms":
                atom_section_start = i + 2  # 跳过标题行和可能的空行
                break
                
        if atom_section_start is not None:
            for i in range(atom_section_start, len(lmp_contents)):
                if lmp_contents[i].strip() == '' or lmp_contents[i].strip().startswith('Bonds'):
                    atom_section_end = i
                    break
        
        # 如果没有找到结束标记，则使用文件末尾
        if atom_section_end is None:
            atom_section_end = len(lmp_contents)
        
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
        atom_section = lmp_contents[atom_section_start:atom_section_end]
        atom_lines = [line for line in atom_section if line.strip() and not line.strip().startswith('#')]
        
        if len(atom_lines) != len(chg_values):
            raise ValueError(f"原子数量({len(atom_lines)})与电荷值数量({len(chg_values)})不匹配")
        
        modified_section = []
        charge_idx = 0
        
        for line in lmp_contents[atom_section_start:atom_section_end]:
            if line.strip() and not line.strip().startswith('#'):
                parts = line.split()
                if len(parts) >= 4:  # 确保行有足够的字段
                    # LAMMPS文件格式：atom-ID atom-type q x y z
                    # 电荷通常是第4个字段（索引3）
                    parts[3] = str(chg_values[charge_idx])
                    modified_line = ' '.join(parts) + '\n'
                    modified_section.append(modified_line)
                    charge_idx += 1
                else:
                    modified_section.append(line)
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
            FileNotFoundError: 当文件不存在时
            ValueError: 当原子数量与电荷值数量不匹配时
        """
        try:
            # 检查文件是否存在
            if not os.path.exists(lmp_file_path):
                raise FileNotFoundError(f"LAMMPS文件不存在: {lmp_file_path}")
            
            if not os.path.exists(chg_file_path):
                raise FileNotFoundError(f"电荷文件不存在: {chg_file_path}")
            
            # 读取文件内容
            lmp_contents = self.read_file(lmp_file_path)
            chg_contents = self.read_file(chg_file_path)
            
            # 提取电荷值
            chg_values = self.extract_charges(chg_contents)
            
            if not chg_values:
                raise ValueError(f"从电荷文件中未能提取到有效的电荷值")
            
            # 查找原子部分位置
            atom_section_start, atom_section_end = self.find_atoms_section(lmp_contents)
            
            if atom_section_start is None or atom_section_end is None:
                raise ValueError(f"在LAMMPS文件中未能找到Atoms部分")
            
            # 修改电荷值
            modified_lmp_contents = self.modify_charges_in_lmp(
                lmp_contents, atom_section_start, atom_section_end, chg_values
            )
            
            # 写入修改后的内容到新文件
            self.write_file(output_lmp_file_path, modified_lmp_contents)
            
            self.logger.info(f"成功修改LAMMPS文件中的电荷值，已保存到: {output_lmp_file_path}")
            
        except Exception as e:
            self.logger.error(f"修改LAMMPS文件中的电荷值时出错: {str(e)}")
            raise 