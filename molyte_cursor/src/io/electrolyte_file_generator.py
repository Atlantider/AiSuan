"""
电解液文件生成器模块，负责处理网页配方数据并生成电解液结构和LAMMPS输入文件
"""
import os
import shutil
import yaml
import subprocess
import time
from pathlib import Path
import tempfile
from ..utils.logger import Logger
from .lammps_generator import LAMMPSGenerator
from .lammps_file_generator import LAMMPSFileGenerator

class ElectrolyteFileGenerator:
    """电解液文件生成器类，用于处理网页配方数据并生成相应的模型结构和输入文件"""
    
    def __init__(self, config_path=None):
        """初始化电解液文件生成器
        
        Args:
            config_path: 配置文件路径，默认使用molyte_cursor/config/paths/default_paths.yaml
        """
        self.logger = Logger().get_logger()
        
        # 设置BOSS环境变量
        boss_dir = "/public/software/boss"
        if os.path.exists(boss_dir):
            os.environ["BOSSdir"] = boss_dir
            self.logger.info(f"已设置BOSSdir环境变量: {boss_dir}")
        else:
            self.logger.warning(f"BOSS软件目录不存在: {boss_dir}")
        
        # 加载配置文件
        if config_path is None:
            # 查找默认配置文件路径
            base_path = Path(__file__).parent.parent.parent
            config_path = base_path / "config" / "paths" / "default_paths.yaml"
        
        self.config = {}
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                self.config = yaml.safe_load(f) or {}
        
        # 获取外部工具路径
        self.paths = self.config.get('paths', {})
        
        # 获取initial_salts目录路径
        self.initial_salts_path = self.paths.get('initial_salts', str(base_path / "inital_salts"))
        if not os.path.exists(self.initial_salts_path):
            os.makedirs(self.initial_salts_path, exist_ok=True)
            self.logger.warning(f"创建initial_salts目录: {self.initial_salts_path}")
        
        # 处理LigParGen路径
        ligpargen_path = self.paths.get('ligpargen')
        if ligpargen_path and os.path.isdir(ligpargen_path):
            self.ligpargen_path = os.path.join(ligpargen_path, "ligpargen")
        else:
            self.ligpargen_path = self.paths.get('ligpargen', 'ligpargen')
        
        self.packmol_path = self.paths.get('packmol', 'packmol')
        self.moltemplate_path = self.paths.get('moltemplate', 'moltemplate.sh')
        
        # 记录找到的工具路径
        self.logger.info(f"使用的外部工具路径: LigParGen={self.ligpargen_path}, Packmol={self.packmol_path}, Moltemplate={self.moltemplate_path}")
        self.logger.info(f"离子文件目录: {self.initial_salts_path}")
        
        # 初始化LAMMPS文件生成器
        self.lammps_generator = LAMMPSGenerator()
    
    def generate_molecule_files(self, molecule_name, smile, output_dir, simulate_mode=False):
        """生成分子文件（包括PDB和力场文件）
        
        Args:
            molecule_name: 分子名称
            smile: SMILES格式的分子结构
            output_dir: 输出目录
            simulate_mode: 是否使用模拟模式（不实际调用外部工具）
            
        Returns:
            生成的文件路径字典
        """
        self.logger.info(f"为分子 {molecule_name} 生成文件 (SMILE: {smile})")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        if simulate_mode:
            # 模拟模式下，只创建空文件
            pdb_path = os.path.join(output_dir, f"{molecule_name}.pdb")
            lt_path = os.path.join(output_dir, f"{molecule_name}.lt")
            
            with open(pdb_path, 'w') as f:
                f.write(f"# 模拟 PDB 文件 for {molecule_name}\n")
            
            with open(lt_path, 'w') as f:
                f.write(f"# 模拟 LT 文件 for {molecule_name}\n")
            
            return {
                'pdb': pdb_path,
                'lt': lt_path
            }
        
        # 真实模式：调用LigParGen生成分子文件
        try:
            # 使用唯一的临时目录名，避免冲突
            temp_dir = os.path.join(output_dir, f"temp_{molecule_name}_{int(time.time())}")
            os.makedirs(temp_dir, exist_ok=True)
            self.logger.info(f"创建临时工作目录: {temp_dir}")
            
            # 准备LigParGen参数
            # 检查是否为离子
            charge = 0
            if "+" in smile:
                charge = 1
            elif "-" in smile:
                charge = -1
            
            # 使用正确的LigParGen参数
            # 根据查看帮助文档，使用 -s 参数指定SMILES，-n 指定分子名称
            cmd = [
                self.ligpargen_path,
                "-s", smile,                   # 指定SMILES字符串
                "-n", molecule_name,           # 指定分子名称
                "-r", "MOL",                   # 设置残基名称
                "-p", temp_dir,                # 输出路径
                "-c", str(charge),             # 设置电荷
                "-o", "0",                     # 优化选项 (0: 不优化，快速)
                "-cgen", "CM1A-LBCC",          # 使用LBCC电荷生成算法
                "-verbose"                     # 输出详细信息
            ]
            
            self.logger.info(f"执行命令: {' '.join(cmd)}")
            process = subprocess.run(cmd, capture_output=True, text=True)
            
            # 输出详细的错误信息和输出信息，帮助调试
            if process.stdout:
                self.logger.info(f"LigParGen输出: {process.stdout}")
            if process.stderr:
                self.logger.warning(f"LigParGen错误: {process.stderr}")
            
            if process.returncode != 0:
                self.logger.error(f"LigParGen执行失败，返回码: {process.returncode}")
                raise Exception(f"LigParGen执行失败: {process.stderr}")
            
            # 列出临时目录中的所有文件，帮助调试
            self.logger.info(f"LigParGen生成的文件:")
            for file_path in Path(temp_dir).glob("*"):
                self.logger.info(f"  - {file_path}")
            
            # 查找生成的文件
            # 检查多种可能的文件命名模式
            possible_pdb_names = [
                f"{molecule_name}.pdb",
                f"{molecule_name}_NEW.pdb",
                f"{molecule_name.upper()}.pdb",
                f"{molecule_name.lower()}.pdb",
                f"MOL.pdb",
                "output.pdb"
            ]
            
            # 检查多种可能的LT文件命名或相关文件
            possible_lt_or_related_names = [
                f"{molecule_name}.lt",
                f"{molecule_name}_NEW.lt",
                f"{molecule_name.upper()}.lt",
                f"{molecule_name.lower()}.lt",
                f"MOL.lt",
                f"{molecule_name}.z",
                f"{molecule_name}_NEW.z",
                f"{molecule_name.upper()}.z",
                f"{molecule_name.lower()}.z",
                f"MOL.z",
                "output.z",
                f"{molecule_name}.xyz",
                f"{molecule_name}_NEW.xyz",
                f"{molecule_name.upper()}.xyz",
                f"{molecule_name.lower()}.xyz"
            ]
            
            # 尝试找到匹配的PDB文件
            pdb_file = None
            for pdb_name in possible_pdb_names:
                file_path = os.path.join(temp_dir, pdb_name)
                if os.path.exists(file_path):
                    pdb_file = file_path
                    self.logger.info(f"找到PDB文件: {pdb_file}")
                    break
            
            # 如果找不到特定命名的文件，尝试找到任何PDB文件
            if not pdb_file:
                pdb_files = list(Path(temp_dir).glob("*.pdb"))
                if pdb_files:
                    pdb_file = str(pdb_files[0])
                    self.logger.info(f"找到PDB文件: {pdb_file}")
                else:
                    # 尝试搜索子目录
                    pdb_files = []
                    for subdir in Path(temp_dir).glob("*"):
                        if subdir.is_dir():
                            pdb_files.extend(list(subdir.glob("*.pdb")))
                    
                    if pdb_files:
                        pdb_file = str(pdb_files[0])
                        self.logger.info(f"在子目录中找到PDB文件: {pdb_file}")
                    else:
                        self.logger.error(f"未找到任何PDB文件")
                        
                        # 创建简化的PDB文件
                        self.logger.warning("创建简化的PDB文件替代LigParGen输出")
                        pdb_file = os.path.join(temp_dir, f"{molecule_name}.pdb")
                        with open(pdb_file, 'w') as f:
                            f.write(f"HETATM    1  C   {molecule_name[:3]} A   1       0.000   0.000   0.000  1.00  0.00           C\n")
                            f.write(f"HETATM    2  O   {molecule_name[:3]} A   1       1.200   0.000   0.000  1.00  0.00           O\n")
                            f.write("END")
                        self.logger.info(f"已创建简化的PDB文件: {pdb_file}")
            
            # 尝试找到匹配的LT或相关文件
            lt_file = None
            z_file = None
            for lt_name in possible_lt_or_related_names:
                file_path = os.path.join(temp_dir, lt_name)
                if os.path.exists(file_path):
                    if lt_name.endswith('.lt'):
                        lt_file = file_path
                        self.logger.info(f"找到LT文件: {lt_file}")
                        break
                    elif lt_name.endswith('.z') and not z_file:
                        z_file = file_path
                        self.logger.info(f"找到Z文件: {z_file}")
            
            # 如果找不到LT文件但找到了Z文件，则从Z文件生成LT文件
            if not lt_file and z_file:
                lt_file = os.path.join(temp_dir, f"{molecule_name}.lt")
                self.logger.info(f"从Z文件生成LT文件: {z_file} -> {lt_file}")
                
                # 读取Z文件内容
                with open(z_file, 'r') as z_f:
                    z_content = z_f.readlines()
                
                # 提取原子和键信息
                atoms_section = False
                bonds_section = False
                atoms_data = []
                bonds_data = []
                
                for line in z_content:
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('ATOMS'):
                        atoms_section = True
                        continue
                    elif line.startswith('BONDS'):
                        atoms_section = False
                        bonds_section = True
                        continue
                    
                    if atoms_section:
                        atoms_data.append(line)
                    elif bonds_section:
                        bonds_data.append(line)
                
                # 生成LT文件内容
                lt_content = [f"# 从Z文件生成的{molecule_name}力场文件"]
                lt_content.append('import "oplsaa.lt"  # 使用OPLS-AA力场')
                lt_content.append('')
                lt_content.append(f'{molecule_name} inherits OPLSAA {{')
                lt_content.append('  # 定义分子结构')
                lt_content.append('  write("Data Atoms") {')
                
                # 添加原子信息
                for i, atom_line in enumerate(atoms_data):
                    parts = atom_line.split()
                    if len(parts) >= 6:
                        atom_id = parts[0]
                        atom_type = parts[1]
                        charge = parts[2]
                        x, y, z = parts[3:6]
                        lt_content.append(f'    $atom:{atom_id} $mol:. @atom:{atom_type} {charge} {x} {y} {z}')
                
                lt_content.append('  }')
                
                # 添加键信息
                if bonds_data:
                    lt_content.append('')
                    lt_content.append('  write("Data Bonds") {')
                    for i, bond_line in enumerate(bonds_data):
                        parts = bond_line.split()
                        if len(parts) >= 3:
                            atom1 = parts[0]
                            atom2 = parts[1]
                            bond_type = parts[2]
                            lt_content.append(f'    $bond:b{i+1} @bond:{bond_type} $atom:{atom1} $atom:{atom2}')
                    lt_content.append('  }')
                
                lt_content.append('}')
                
                # 写入LT文件
                with open(lt_file, 'w') as lt_f:
                    lt_f.write('\n'.join(lt_content))
            
            # 如果仍然没有LT文件，创建一个简化的LT文件
            if not lt_file:
                lt_file = os.path.join(temp_dir, f"{molecule_name}.lt")
                self.logger.warning(f"未找到任何LT或Z文件，创建简化的LT文件")
                
                lt_content = f"""# 简化的{molecule_name}力场文件
import "oplsaa.lt"  # 使用OPLS-AA力场

{molecule_name} inherits OPLSAA {{
  # 定义简化分子结构（基于基本的醚结构）
  write("Data Atoms") {{
    $atom:C1 $mol:. @atom:90 0.0 0.0 0.0 0.0  # 碳原子
    $atom:O1 $mol:. @atom:180 0.0 1.0 0.0 0.0  # 氧原子
    $atom:C2 $mol:. @atom:90 0.0 2.0 0.0 0.0  # 碳原子
  }}
  
  write("Data Bonds") {{
    $bond:b1 @bond:30 $atom:C1 $atom:O1  # C-O键
    $bond:b2 @bond:30 $atom:O1 $atom:C2  # O-C键
  }}
}}
"""
                with open(lt_file, 'w') as f:
                    f.write(lt_content)
                self.logger.info(f"已创建简化的LT文件: {lt_file}")
            
            # 复制文件到输出目录
            pdb_path = os.path.join(output_dir, f"{molecule_name}.pdb")
            lt_path = os.path.join(output_dir, f"{molecule_name}.lt")
            
            shutil.copy2(pdb_file, pdb_path)
            shutil.copy2(lt_file, lt_path)
            
            self.logger.info(f"已复制文件到输出目录: {pdb_path}, {lt_path}")
            
            # 不删除临时目录，以便稍后检查（在生产环境可以打开清理）
            # shutil.rmtree(temp_dir)
            
            return {
                'pdb': pdb_path,
                'lt': lt_path
            }
        except Exception as e:
            self.logger.error(f"生成分子文件时出错: {str(e)}")
            # 出错时创建简化文件
            self.logger.warning("出错时创建简化文件")
            pdb_path = os.path.join(output_dir, f"{molecule_name}.pdb")
            lt_path = os.path.join(output_dir, f"{molecule_name}.lt")
            
            with open(pdb_path, 'w') as f:
                f.write(f"HETATM    1  C   {molecule_name[:3]} A   1       0.000   0.000   0.000  1.00  0.00           C\n")
                f.write(f"HETATM    2  O   {molecule_name[:3]} A   1       1.200   0.000   0.000  1.00  0.00           O\n")
                f.write("END")
            
            with open(lt_path, 'w') as f:
                f.write(f"""# 简化的{molecule_name}力场文件
import "oplsaa.lt"  # 使用OPLS-AA力场

{molecule_name} inherits OPLSAA {{
  # 定义简化分子结构（基于基本的醚结构）
  write("Data Atoms") {{
    $atom:C1 $mol:. @atom:90 0.0 0.0 0.0 0.0  # 碳原子
    $atom:O1 $mol:. @atom:180 0.0 1.0 0.0 0.0  # 氧原子
    $atom:C2 $mol:. @atom:90 0.0 2.0 0.0 0.0  # 碳原子
  }}
  
  write("Data Bonds") {{
    $bond:b1 @bond:30 $atom:C1 $atom:O1  # C-O键
    $bond:b2 @bond:30 $atom:O1 $atom:C2  # O-C键
  }}
}}
""")
            
            return {
                'pdb': pdb_path,
                'lt': lt_path
            }
    
    def run_packmol(self, system_name, molecule_files, box_size, output_dir):
        """运行Packmol生成初始结构
        
        Args:
            system_name: 系统名称
            molecule_files: 分子文件路径字典（名称到PDB文件路径的映射）
            box_size: 模拟盒子大小（Å）
            output_dir: 输出目录
            
        Returns:
            生成的PDB文件路径
            
        Raises:
            FileNotFoundError: 当找不到Packmol可执行文件或输入文件时
            RuntimeError: 当Packmol执行失败时
        """
        self.logger.info(f"为系统 {system_name} 运行Packmol")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 准备Packmol输入文件
        inp_path = os.path.join(output_dir, f"{system_name}.inp")
        with open(inp_path, 'w') as f:
            f.write(f"tolerance 2.0\nfiletype pdb\noutput {system_name}.pdb\n\n")
            
            # 为每个分子添加结构信息
            for molecule_name, molecule_info in molecule_files.items():
                number = molecule_info.get('number', 1)
                pdb_path = molecule_info.get('pdb', '')
                
                if not pdb_path or not os.path.exists(pdb_path):
                    error_msg = f"未找到分子 {molecule_name} 的PDB文件: {pdb_path}"
                    self.logger.error(error_msg)
                    raise FileNotFoundError(error_msg)
                
                # 使用绝对路径确保Packmol可以找到文件
                pdb_path_abs = os.path.abspath(pdb_path)
                f.write(f"structure {pdb_path_abs}\n")
                f.write(f"  number {number}\n")
                f.write(f"  inside box 0.0 0.0 0.0 {box_size} {box_size} {box_size}\n")
                f.write("end structure\n\n")
        
        output_pdb = os.path.join(output_dir, f"{system_name}.pdb")
        
        # 检查Packmol路径
        packmol_path = self.packmol_path
        if not os.path.exists(packmol_path) and not os.path.isabs(packmol_path):
            # 如果提供的不是绝对路径且文件不存在，尝试查找完整路径
            packmol_path = shutil.which("packmol") or "/public/software/packmol-20.16.0/packmol"
            
            if not os.path.exists(packmol_path):
                error_msg = f"找不到packmol可执行文件，请检查配置或安装"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
        
        # 切换到输出目录，因为Packmol会在当前目录生成输出文件
        current_dir = os.getcwd()
        os.chdir(output_dir)
        
        try:
            # 使用直接的命令行参数而不是输入重定向
            cmd = [packmol_path, "-i", os.path.basename(inp_path)]
            
            self.logger.info(f"执行命令: {' '.join(cmd)}")
            self.logger.info(f"工作目录: {output_dir}")
            
            # 不使用shell模式执行命令
            process = subprocess.run(cmd, capture_output=True, text=True)
            
            # 检查命令是否成功
            if process.returncode != 0:
                # 如果命令失败，尝试另一种方法
                self.logger.warning(f"标准Packmol命令失败，尝试备用方法: {process.stderr}")
                
                # 尝试备用方法：直接运行可执行文件并将输入文件读入内存传递
                cmd_alt = [packmol_path]
                
                with open(inp_path, 'r') as f:
                    inp_content = f.read()
                
                process = subprocess.run(cmd_alt, input=inp_content, capture_output=True, text=True)
                
                if process.returncode != 0:
                    error_msg = f"Packmol执行失败: {process.stderr}"
                    self.logger.error(error_msg)
                    raise RuntimeError(error_msg)
            
            # 检查输出文件是否存在
            if not os.path.exists(os.path.join(output_dir, f"{system_name}.pdb")):
                error_msg = f"Packmol未生成输出文件: {output_pdb}"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
            
            # 输出Packmol的执行结果
            if process.stdout:
                self.logger.info(f"Packmol输出: {process.stdout}")
            
            self.logger.info(f"Packmol成功生成PDB文件: {output_pdb}")
            return output_pdb
        
        finally:
            # 恢复原来的工作目录（即使出现异常也会执行）
            os.chdir(current_dir)
    
    def run_moltemplate(self, system_name, lt_files, box_size, output_dir):
        """运行Moltemplate生成LAMMPS数据文件
        
        Args:
            system_name: 系统名称
            lt_files: LT文件路径字典
            box_size: 模拟盒子大小（Å）
            output_dir: 输出目录
            
        Returns:
            生成的LAMMPS数据文件路径
            
        Raises:
            FileNotFoundError: 当找不到Moltemplate可执行文件或输入文件时
            RuntimeError: 当Moltemplate执行失败时
        """
        self.logger.info(f"为系统 {system_name} 运行Moltemplate")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 准备系统LT文件
        system_lt_path = os.path.join(output_dir, f"{system_name}.lt")
        with open(system_lt_path, 'w') as f:
            # 导入所有分子LT文件
            for molecule_name, lt_path in lt_files.items():
                if isinstance(lt_path, dict):
                    lt_path = lt_path.get('lt', '')
                
                if not lt_path or not os.path.exists(lt_path):
                    error_msg = f"未找到分子 {molecule_name} 的LT文件: {lt_path}"
                    self.logger.error(error_msg)
                    raise FileNotFoundError(error_msg)
                
                f.write(f'import "{lt_path}"\n')
            
            f.write("\n")
            
            # 实例化分子
            for molecule_name, lt_info in lt_files.items():
                if isinstance(lt_info, dict):
                    number = lt_info.get('number', 1)
                else:
                    number = 1
                f.write(f'{molecule_name}s = new {molecule_name}[{number}]\n')
            
            f.write("\n")
            
            # 添加模拟盒子边界定义
            f.write(f"""write_once("Data Boundary") {{
0 {box_size} xlo xhi
0 {box_size} ylo yhi
0 {box_size} zlo zhi
}}
""")
        
        output_data = os.path.join(output_dir, f"{system_name}.data")
        
        # 检查Moltemplate路径
        moltemplate_path = self.moltemplate_path
        if not os.path.exists(moltemplate_path) and not os.path.isabs(moltemplate_path):
            # 如果提供的不是绝对路径且文件不存在，尝试查找完整路径
            moltemplate_path = shutil.which("moltemplate.sh") or "/public/software/moltemplate/moltemplate.sh"
            
            if not os.path.exists(moltemplate_path):
                error_msg = f"找不到moltemplate可执行文件，请检查配置或安装"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
        
        # 验证PDB文件是否存在
        pdb_path = os.path.join(output_dir, f"{system_name}.pdb")
        if not os.path.exists(pdb_path):
            error_msg = f"未找到PDB文件: {pdb_path}，Moltemplate需要该文件"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        # 保存当前目录
        current_dir = os.getcwd()
        
        try:
            # 切换到输出目录
            os.chdir(output_dir)
            
            # 构建命令
            cmd = [
                moltemplate_path, 
                "-xyz", os.path.basename(pdb_path), 
                os.path.basename(system_lt_path)
            ]
            
            self.logger.info(f"执行命令: {' '.join(cmd)}")
            self.logger.info(f"工作目录: {output_dir}")
            
            # 执行Moltemplate命令
            process = subprocess.run(cmd, capture_output=True, text=True)
            
            # 输出Moltemplate的结果
            if process.stdout:
                self.logger.info(f"Moltemplate输出: {process.stdout}")
            if process.stderr:
                self.logger.warning(f"Moltemplate错误: {process.stderr}")
            
            if process.returncode != 0:
                error_msg = f"Moltemplate执行失败: {process.stderr}"
                self.logger.error(error_msg)
                raise RuntimeError(error_msg)
            
            # 检查生成的文件
            if not os.path.exists(output_data):
                # 尝试找到Moltemplate可能生成的其他数据文件
                data_files = list(Path(output_dir).glob("*.data"))
                if data_files:
                    alternative_data = str(data_files[0])
                    self.logger.info(f"找到替代数据文件: {alternative_data}")
                    shutil.copy2(alternative_data, output_data)
                else:
                    error_msg = f"Moltemplate未生成任何LAMMPS数据文件"
                    self.logger.error(error_msg)
                    raise FileNotFoundError(error_msg)
            
            self.logger.info(f"Moltemplate成功生成数据文件: {output_data}")
            return output_data
        
        finally:
            # 恢复原来的工作目录（即使出现异常也会执行）
            os.chdir(current_dir)
    
    def generate_inp_file(self, system_name, molecule_components, box_size, output_dir):
        """生成Packmol输入文件
        
        Args:
            system_name: 系统名称
            molecule_components: 分子组分列表
            box_size: 模拟盒子大小
            output_dir: 输出目录
            
        Returns:
            生成的INP文件路径
        """
        self.logger.info(f"为系统 {system_name} 生成INP文件")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 计算每种组分的分子数量
        # 假设密度为1.0 g/cm³，使用摩尔质量和摩尔浓度计算
        # 这里采用简化计算方法
        box_volume = box_size ** 3  # Å³
        box_volume_L = box_volume * 1e-27  # 转换为L
        
        # 准备分子文件路径字典
        molecule_files = {}
        
        # 创建inp文件
        inp_file_path = os.path.join(output_dir, f"{system_name}.inp")
        with open(inp_file_path, 'w') as f:
            f.write(f"# Packmol input file for {system_name}\n")
            f.write("tolerance 2.0\n")
            f.write(f"output {system_name}.pdb\n")
            f.write(f"filetype pdb\n\n")
            
            for component in molecule_components:
                name = component['name']
                number = component['number']
                pdb_file = component.get('pdb_file')
                
                if not pdb_file or not os.path.exists(pdb_file):
                    self.logger.warning(f"未找到PDB文件: {pdb_file}")
                    continue
                
                # 使用绝对路径，确保Packmol可以找到文件
                pdb_file_abs = os.path.abspath(pdb_file)
                
                f.write(f"# Component: {name}, Number: {number}\n")
                f.write(f"structure {pdb_file_abs}\n")
                f.write(f"  number {number}\n")
                f.write(f"  inside box 0.0 0.0 0.0 {box_size} {box_size} {box_size}\n")
                f.write("end structure\n\n")
                
                # 添加到分子文件字典
                molecule_files[name] = {
                    'number': number,
                    'pdb': pdb_file_abs
                }
        
        self.logger.info(f"已生成INP文件: {inp_file_path}")
        return inp_file_path, molecule_files

    def run_ligpargen(self, smile, charge, molecule_name, output_dir):
        """运行LigParGen生成分子结构文件
        
        Args:
            smile: 分子的SMILES字符串
            charge: 分子的电荷
            molecule_name: 分子名称
            output_dir: 输出目录
            
        Returns:
            文件路径字典: {'.pdb': pdb路径, '.lt': lt路径}
            
        Raises:
            FileNotFoundError: 当LigParGen未找到或执行失败时
            RuntimeError: 当转换过程出错时
        """
        self.logger.info(f"为分子 {molecule_name} 运行LigParGen")
        
        try:
            # 创建临时目录
            temp_dir = tempfile.mkdtemp()
            
            # 设置LigParGen命令
            cmd = [
                self.ligpargen_path,
                "-s", smile,            # SMILES字符串
                "-n", molecule_name,    # 分子名称
                "-r", "MOL",            # 输出格式
                "-p", temp_dir,         # 输出路径
                "-c", str(charge),      # 分子电荷
                "-o", "0",              # 优化级别 (0-3)
                "-cgen", "CM1A-LBCC",   # 电荷生成方法
                "-verbose"              # 详细输出
            ]
            
            # 执行LigParGen命令
            self.logger.info(f"执行LigParGen命令: {' '.join(cmd)}")
            process = subprocess.run(cmd, check=False, capture_output=True, text=True)
            
            if process.returncode == 0:
                self.logger.info("LigParGen命令执行成功")
                self.logger.info(f"LigParGen输出: {process.stdout}")
            else:
                self.logger.error(f"LigParGen执行失败，退出码: {process.returncode}")
                self.logger.error(f"错误信息: {process.stderr}")
                raise RuntimeError(f"LigParGen执行失败: {process.stderr}")
            
            # 列出临时目录中的文件
            self.logger.info(f"检查临时目录中的文件: {temp_dir}")
            found_files = os.listdir(temp_dir)
            
            for file in found_files:
                self.logger.info(f"找到文件: {file}")
            
            # 获取PDB文件
            pdb_file = os.path.join(temp_dir, f"{molecule_name}.openmm.pdb")
            if not os.path.exists(pdb_file):
                # 尝试其他可能的PDB文件
                alternative_pdb = os.path.join(temp_dir, f"{molecule_name}.q.pdb")
                if os.path.exists(alternative_pdb):
                    pdb_file = alternative_pdb
                else:
                    self.logger.error(f"找不到PDB文件: {pdb_file} 或 {alternative_pdb}")
                    raise FileNotFoundError(f"LigParGen未生成PDB文件")
            
            # 查找LAMMPS文件
            lammps_file = os.path.join(temp_dir, f"{molecule_name}.lammps.lmp")
            if not os.path.exists(lammps_file):
                self.logger.error(f"找不到LAMMPS文件: {lammps_file}")
                raise FileNotFoundError(f"LigParGen未生成LAMMPS文件")
            
            self.logger.info(f"找到LAMMPS文件: {lammps_file}")
            
            # 创建输出目录
            os.makedirs(output_dir, exist_ok=True)
            
            # 目标文件路径
            dest_pdb = os.path.join(output_dir, f"{molecule_name}.pdb")
            dest_lt = os.path.join(output_dir, f"{molecule_name}.lt")
            
            # 复制PDB文件
            shutil.copy2(pdb_file, dest_pdb)
            
            # 将LAMMPS文件转换为LT文件
            self.convert_lammps_to_lt(lammps_file, dest_lt, molecule_name)
            
            self.logger.info(f"已复制文件到输出目录: {dest_pdb}")
            self.logger.info(f"已生成LT文件: {dest_lt}")
            
            # 清理临时目录
            shutil.rmtree(temp_dir)
            
            return {'.pdb': dest_pdb, '.lt': dest_lt}
            
        except Exception as e:
            self.logger.error(f"LigParGen处理过程出错: {str(e)}")
            raise FileNotFoundError(f"LigParGen处理失败: {str(e)}")
    
    def convert_lammps_to_lt(self, lammps_file, lt_file, molecule_name):
        """将LAMMPS数据文件转换为LT文件
        
        Args:
            lammps_file: LAMMPS数据文件路径
            lt_file: 输出LT文件路径
            molecule_name: 分子名称
        
        Raises:
            RuntimeError: 当转换过程出错时
        """
        self.logger.info(f"将LAMMPS文件转换为LT文件: {lammps_file} -> {lt_file}")
        
        try:
            # 解析LAMMPS文件
            with open(lammps_file, 'r') as f:
                lammps_content = f.readlines()
            
            # 提取原子、键、角度、二面角信息
            atoms_section = False
            bonds_section = False
            angles_section = False
            dihedrals_section = False
            
            atoms = []
            bonds = []
            angles = []
            dihedrals = []
            
            # 检测数据部分
            for line in lammps_content:
                line = line.strip()
                
                # 检测部分开始
                if "Atoms" in line:
                    atoms_section = True
                    continue
                elif "Bonds" in line:
                    atoms_section = False
                    bonds_section = True
                    continue
                elif "Angles" in line:
                    bonds_section = False
                    angles_section = True
                    continue
                elif "Dihedrals" in line:
                    angles_section = False
                    dihedrals_section = True
                    continue
                
                # 空行或注释行跳过
                if not line or line.startswith('#'):
                    continue
                
                # 根据当前部分读取数据
                if atoms_section:
                    parts = line.split()
                    if len(parts) >= 7:  # 原子ID 分子ID 原子类型 电荷 x y z
                        atoms.append(parts)
                elif bonds_section:
                    parts = line.split()
                    if len(parts) >= 4:  # 键ID 键类型 原子1 原子2
                        bonds.append(parts)
                elif angles_section:
                    parts = line.split()
                    if len(parts) >= 5:  # 角ID 角类型 原子1 原子2 原子3
                        angles.append(parts)
                elif dihedrals_section:
                    parts = line.split()
                    if len(parts) >= 6:  # 二面角ID 二面角类型 原子1 原子2 原子3 原子4
                        dihedrals.append(parts)
            
            # 生成LT文件
            with open(lt_file, 'w') as f:
                f.write(f"# 由LAMMPS格式转换的LT文件 - {molecule_name}\n")
                f.write("import \"oplsaa.lt\"\n\n")
                f.write(f"{molecule_name} inherits OPLSAA {{\n")
                
                # 写入原子部分
                f.write("  # 原子\n")
                f.write("  write(\"Data Atoms\") {\n")
                for atom in atoms:
                    if len(atom) >= 7:
                        atom_id = atom[0]
                        # mol_id = atom[1]  # 分子ID在LT中通常不需要
                        atom_type = atom[2]
                        charge = atom[3]
                        x, y, z = atom[4], atom[5], atom[6]
                        f.write(f"    $atom:{atom_id} $mol:. @atom:{atom_type} {charge} {x} {y} {z}\n")
                f.write("  }\n\n")
                
                # 写入键部分
                if bonds:
                    f.write("  # 键\n")
                    f.write("  write(\"Data Bonds\") {\n")
                    for bond in bonds:
                        if len(bond) >= 4:
                            bond_id = bond[0]
                            bond_type = bond[1]
                            atom1, atom2 = bond[2], bond[3]
                            f.write(f"    $bond:{bond_id} @bond:{bond_type} $atom:{atom1} $atom:{atom2}\n")
                    f.write("  }\n\n")
                
                # 写入角部分
                if angles:
                    f.write("  # 角\n")
                    f.write("  write(\"Data Angles\") {\n")
                    for angle in angles:
                        if len(angle) >= 5:
                            angle_id = angle[0]
                            angle_type = angle[1]
                            atom1, atom2, atom3 = angle[2], angle[3], angle[4]
                            f.write(f"    $angle:{angle_id} @angle:{angle_type} $atom:{atom1} $atom:{atom2} $atom:{atom3}\n")
                    f.write("  }\n\n")
                
                # 写入二面角部分
                if dihedrals:
                    f.write("  # 二面角\n")
                    f.write("  write(\"Data Dihedrals\") {\n")
                    for dihedral in dihedrals:
                        if len(dihedral) >= 6:
                            dih_id = dihedral[0]
                            dih_type = dihedral[1]
                            atom1, atom2, atom3, atom4 = dihedral[2], dihedral[3], dihedral[4], dihedral[5]
                            f.write(f"    $dihedral:{dih_id} @dihedral:{dih_type} $atom:{atom1} $atom:{atom2} $atom:{atom3} $atom:{atom4}\n")
                    f.write("  }\n\n")
                
                # 结束分子定义
                f.write("}\n")
            
            self.logger.info(f"成功将LAMMPS文件转换为LT文件: {lt_file}")
            
        except Exception as e:
            error_msg = f"LAMMPS到LT转换失败: {str(e)}"
            self.logger.error(error_msg)
            raise RuntimeError(error_msg)

    def create_ion_files(self, ion_name, charge, output_dir):
        """为离子创建简化的PDB和LT文件
        
        Args:
            ion_name: 离子名称
            charge: 离子电荷
            output_dir: 输出目录
            
        Returns:
            文件路径字典: {'.pdb': pdb路径, '.lt': lt路径}
        """
        self.logger.info(f"为离子 {ion_name} 创建简化结构文件 (电荷: {charge})")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 处理离子名称，移除+或-符号
        clean_ion_name = ion_name.replace('+', '').replace('-', '')
        
        # 准备输出文件路径
        pdb_file = os.path.join(output_dir, f"{clean_ion_name}.pdb")
        lt_file = os.path.join(output_dir, f"{clean_ion_name}.lt")
        
        # 检查文件是否已存在
        if os.path.exists(pdb_file) and os.path.exists(lt_file):
            self.logger.info(f"离子文件已存在，跳过创建: {pdb_file}, {lt_file}")
            return {'.pdb': pdb_file, '.lt': lt_file}
        
        # 创建简化的PDB文件
        with open(pdb_file, 'w') as f:
            f.write(f"TITLE     简化的 {ion_name} 离子结构\n")
            f.write(f"REMARK    电荷: {charge}\n")
            f.write(f"HETATM    1  {clean_ion_name[:1]:1s}   {clean_ion_name:3s} A   1       0.000   0.000   0.000  1.00  0.00           {clean_ion_name[:1]:1s}\n")
            f.write("END\n")
        
        # 创建简化的LT文件
        with open(lt_file, 'w') as f:
            f.write(f"""# 简化的 {ion_name} 离子力场文件
import "oplsaa.lt"  # 使用OPLS-AA力场

{clean_ion_name} inherits OPLSAA {{
  # 定义简化的离子结构
  write("Data Atoms") {{
    $atom:1 $mol:. @atom:{406 if charge > 0 else 407} {charge:.1f} 0.0 0.0 0.0  # 离子中心原子
  }}
}}
""")
        
        self.logger.info(f"成功创建离子文件: {pdb_file}, {lt_file}")
        return {'.pdb': pdb_file, '.lt': lt_file}

    def load_ion_files(self, ion_name, output_dir):
        """从initial_salts目录加载离子文件
        
        Args:
            ion_name: 离子名称
            output_dir: 输出目录
            
        Returns:
            文件路径字典: {'.pdb': pdb路径, '.lt': lt路径}
            
        Raises:
            FileNotFoundError: 当找不到离子文件时
        """
        self.logger.info(f"从initial_salts目录加载离子 {ion_name} 的文件")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 处理离子名称，移除+或-符号
        clean_ion_name = ion_name.replace('+', '').replace('-', '')
        
        # 在initial_salts目录中查找文件
        source_pdb = os.path.join(self.initial_salts_path, f"{clean_ion_name}.pdb")
        source_lt = os.path.join(self.initial_salts_path, f"{clean_ion_name}.lt")
        
        # 检查文件是否存在
        if not os.path.exists(source_pdb) or not os.path.exists(source_lt):
            error_msg = f"在initial_salts目录中找不到离子文件: {source_pdb} 或 {source_lt}"
            self.logger.error(error_msg)
            
            # 如果找不到文件，尝试创建简化的文件作为备选
            self.logger.warning(f"找不到离子文件，将使用简化模型")
            return self.create_ion_files(clean_ion_name, 1 if clean_ion_name in ["Li", "Na", "K"] else -1, output_dir)
        
        # 目标文件路径，也使用不带电荷符号的文件名
        dest_pdb = os.path.join(output_dir, f"{clean_ion_name}.pdb")
        dest_lt = os.path.join(output_dir, f"{clean_ion_name}.lt")
        
        # 复制文件
        shutil.copy2(source_pdb, dest_pdb)
        shutil.copy2(source_lt, dest_lt)
        
        self.logger.info(f"成功复制离子文件: {dest_pdb}, {dest_lt}")
        return {'.pdb': dest_pdb, '.lt': dest_lt}

    def generate_lammps_input_files(self, molecule_files, output_dir, packmol_output, config):
        """生成LAMMPS输入文件
        
        将Packmol输出和力场文件合并生成LAMMPS输入文件
        
        Args:
            molecule_files: 分子文件信息
            output_dir: 输出目录
            packmol_output: Packmol输出信息
            config: 配置数据
            
        Returns:
            Dict: 包含生成的文件信息
        """
        self.logger.info("生成LAMMPS输入文件...")
        
        try:
            # 确保输出目录存在
            os.makedirs(output_dir, exist_ok=True)
            
            # 使用LAMMPSFileGenerator生成LAMMPS输入文件
            file_generator = LAMMPSFileGenerator()
            
            # 将Packmol输出和力场文件合并生成LAMMPS输入文件
            lammps_input_path = os.path.join(output_dir, 'lammps.input')
            data_file_path = os.path.join(output_dir, 'system.data')
            
            # 生成LAMMPS输入文件时提供Packmol输出和力场信息
            input_file_path = file_generator.generate_input_file(
                config=config, 
                output_path=lammps_input_path,
                data_file=packmol_output.get('data_file'),
                topology_files=molecule_files.get('topology_files', []),
                forcefield_files=molecule_files.get('forcefield_files', [])
            )
            
            # 创建生成的文件列表
            generated_files = [{
                'type': 'lammps_input',
                'path': lammps_input_path,
                'is_main': True
            }]
            
            # 添加其他生成的文件
            if os.path.exists(data_file_path):
                generated_files.append({
                    'type': 'lammps_data',
                    'path': data_file_path,
                    'is_main': False
                })
            
            # 添加力场文件
            for ff_file in molecule_files.get('forcefield_files', []):
                generated_files.append({
                    'type': 'forcefield',
                    'path': ff_file,
                    'is_main': False
                })
            
            return generated_files
            
        except Exception as e:
            self.logger.error(f"生成LAMMPS输入文件失败: {str(e)}")
            raise

def generate_electrolyte_input_files(recipe, output_dir):
    """从电解液配方生成所有必需的输入文件
    
    Args:
        recipe: ElectrolyteRecipe对象，包含配方的完整信息
        output_dir: 输出目录路径
        
    Returns:
        生成的文件路径字典
    """
    generator = ElectrolyteFileGenerator()
    logger = generator.logger
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 处理所有分子组分
    molecule_files = {}
    
    # 处理阳离子
    for cation in recipe.cations:
        if cation.smile:  # 如果有SMILE，使用LigParGen生成文件
            files = generator.generate_molecule_files(
                molecule_name=cation.name,
                smile=cation.smile,
                output_dir=output_dir
            )
        else:  # 否则从initial_salts目录加载离子文件
            files = generator.load_ion_files(cation.name, output_dir)
        molecule_files[cation.name] = files
    
    # 处理阴离子
    for anion in recipe.anions:
        if anion.smile:
            files = generator.generate_molecule_files(
                molecule_name=anion.name,
                smile=anion.smile,
                output_dir=output_dir
            )
        else:
            files = generator.load_ion_files(anion.name, output_dir)
        molecule_files[anion.name] = files
    
    # 处理溶剂分子
    for solvent in recipe.solvents:
        files = generator.generate_molecule_files(
            molecule_name=solvent.name,
            smile=solvent.smile,
            output_dir=output_dir
        )
        molecule_files[solvent.name] = files
    
    # 生成Packmol输入文件并运行
    packmol_output = generator.run_packmol(
        system_name=recipe.name,
        molecule_files={
            **{cat.name: {'pdb': molecule_files[cat.name]['pdb'], 'number': cat.number}
               for cat in recipe.cations},
            **{an.name: {'pdb': molecule_files[an.name]['pdb'], 'number': an.number}
               for an in recipe.anions},
            **{sol.name: {'pdb': molecule_files[sol.name]['pdb'], 'number': sol.number}
               for sol in recipe.solvents}
        },
        box_size=recipe.box_size,
        output_dir=output_dir
    )
    
    # 生成Moltemplate输入文件并运行
    moltemplate_output = generator.run_moltemplate(
        system_name=recipe.name,
        lt_files={
            **{cat.name: molecule_files[cat.name]['lt'] for cat in recipe.cations},
            **{an.name: molecule_files[an.name]['lt'] for an in recipe.anions},
            **{sol.name: molecule_files[sol.name]['lt'] for sol in recipe.solvents}
        },
        box_size=recipe.box_size,
        output_dir=output_dir
    )
    
    # 生成LAMMPS输入文件
    lammps_output = generator.lammps_generator.generate_lammps_input(
        system_name=recipe.name,
        temperature=recipe.temperature,
        components={
            'cations': [{'name': cat.name, 'number': cat.number} for cat in recipe.cations],
            'anions': [{'name': an.name, 'number': an.number} for an in recipe.anions],
            'solvents': [{'name': sol.name, 'number': sol.number} for sol in recipe.solvents]
        },
        output_dir=output_dir
    )
    
    # 返回所有生成的文件路径
    return {
        'molecule_files': molecule_files,
        'packmol_output': packmol_output,
        'moltemplate_output': moltemplate_output,
        'lammps_output': lammps_output
    } 