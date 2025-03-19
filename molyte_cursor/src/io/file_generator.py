"""
文件生成器模块，负责生成各种配置和输入文件
"""
from pathlib import Path
from ..utils.logger import Logger
import os
import yaml

class FileGenerator:
    """文件生成器类"""
    
    def __init__(self):
        """初始化文件生成器"""
        self.logger = Logger().get_logger()
    
    def create_job_script(self, destination_path, system_name):
        """创建作业脚本
        
        Args:
            destination_path: 目标目录路径
            system_name: 系统名称
            
        Returns:
            作业脚本文件路径
        """
        script_content = f"""#!/bin/bash
#PBS -N {system_name}
#PBS -o out.dat
#PBS -e err.dat
#PBS -q batch
#PBS -l nodes=1:ppn=64
#PBS -l walltime=7200:00:00
cd $PBS_O_WORKDIR

source /home/iei/share/software/environment/lammps.env
EXEC=/home/iei/share/software/lammps/200303/lmp_mpi
input=all.in
mpirun -n 64 $EXEC <{system_name}.in >{system_name}.log
"""
        
        job_script_path = destination_path / 'job.sh'
        with open(job_script_path, 'w') as job_file:
            job_file.write(script_content)
            
        self.logger.info(f"Generated job script: {job_script_path}")
        return job_script_path
    
    def create_packmol_input(self, system, destination_path):
        """创建Packmol输入文件
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            
        Returns:
            Packmol输入文件路径
        """
        packmol_input = f"tolerance 2.0\nfiletype pdb\noutput {system.name}.pdb\n"
        
        # 添加阳离子
        for cation in system.cations:
            if cation.is_valid():
                packmol_input += f"""structure {cation.name}.pdb
   number {cation.number}
   inside box 0. 0. 0 {system.box_size} {system.box_size} {system.box_size}
end structure
"""
        
        # 添加阴离子
        for anion in system.anions:
            if anion.is_valid():
                packmol_input += f"""structure {anion.name}.pdb
   number {anion.number}
   inside box 0. 0. 0 {system.box_size} {system.box_size} {system.box_size}
end structure
"""
        
        # 添加溶剂
        for solvent in system.solvents:
            if solvent.is_valid():
                packmol_input += f"""structure {solvent.name}.charmm.pdb
   number {solvent.number}
   inside box 0. 0. 0 {system.box_size} {system.box_size} {system.box_size}
end structure
"""
        
        packmol_input_path = destination_path / f"{system.name}.inp"
        with open(packmol_input_path, 'w') as packmol_file:
            packmol_file.write(packmol_input)
            
        self.logger.info(f"Generated Packmol input file: {packmol_input_path}")
        return packmol_input_path
    
    def create_moltemplate_lt(self, system, destination_path):
        """创建Moltemplate LT文件
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            
        Returns:
            Moltemplate LT文件路径
        """
        moltemplate_content = ""
        
        # 导入分子LT文件
        for cation in system.cations:
            if cation.is_valid():
                moltemplate_content += f'import "{cation.name}.lt"\n'
                
        for anion in system.anions:
            if anion.is_valid():
                moltemplate_content += f'import "{anion.name}.lt"\n'
                
        for solvent in system.solvents:
            if solvent.is_valid():
                moltemplate_content += f'import "{solvent.name}.lt"\n'
        
        # 实例化分子
        for cation in system.cations:
            if cation.is_valid():
                moltemplate_content += f'{cation.name}s = new {cation.name}[{cation.number}]\n'
                
        for anion in system.anions:
            if anion.is_valid():
                moltemplate_content += f'{anion.name}s = new {anion.name}[{anion.number}]\n'
                
        for solvent in system.solvents:
            if solvent.is_valid():
                moltemplate_content += f'{solvent.name}s = new {solvent.name}[{solvent.number}]\n'
        
        # 添加模拟盒子边界定义
        moltemplate_content += f"""
write_once("Data Boundary") {{
    0 {system.box_size} xlo xhi
    0 {system.box_size} ylo yhi
    0 {system.box_size} zlo zhi
}}
"""
        
        moltemplate_lt_path = destination_path / f"{system.name}.lt"
        with open(moltemplate_lt_path, 'w') as moltemplate_file:
            moltemplate_file.write(moltemplate_content)
            
        self.logger.info(f"Generated Moltemplate LT file: {moltemplate_lt_path}")
        return moltemplate_lt_path
    
    def create_element_list(self, lmp_data_path):
        """从LAMMPS数据文件创建元素列表
        
        Args:
            lmp_data_path: LAMMPS数据文件路径
            
        Returns:
            (元素列表, RDF对列表, RDF标签列表)
        """
        def element_type(atomic_mass):
            """根据原子质量返回元素符号"""
            element_map = {
                1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C',
                14: 'N', 16: 'O', 19: 'F', 20: 'Ne', 23: 'Na', 24: 'Mg',
                27: 'Al', 28: 'Si', 31: 'P', 32: 'S', 35: 'Cl', 39: 'K',
                40: 'Ca', 45: 'Sc', 48: 'Ti', 51: 'V', 52: 'Cr', 55: 'Mn',
                56: 'Fe', 59: 'Co', 64: 'Cu', 65: 'Zn', 70: 'Ga', 73: 'Ge',
                75: 'As', 79: 'Se'
            }
            return element_map.get(atomic_mass, 'A')
        
        # 解析LAMMPS数据文件
        element_list = []
        
        with open(lmp_data_path, 'r', encoding='utf-8') as lmp_file:
            content_list = lmp_file.readlines()
        
        # 查找Masses部分
        for i, line in enumerate(content_list):
            if line.strip().lower() == 'masses':
                for j in range(i + 2, len(content_list)):
                    if not content_list[j].strip() or 'atoms' in content_list[j].lower():
                        break
                    
                    parts = content_list[j].split()
                    if len(parts) < 2:
                        self.logger.warning(f"Malformed Masses line: '{content_list[j]}'")
                        continue
                    
                    try:
                        atomic_mass = round(float(parts[1]))
                        element = element_type(atomic_mass)
                        element_list.append(element)
                    except ValueError:
                        self.logger.warning(f"Invalid atomic mass '{parts[1]}' in Masses line.")
                break
        
        # 生成RDF对和标签
        rdf_pairs = ''
        rdf_labels = ''
        
        # 阳离子元素列表
        cation_elements = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Zn', 'Al']
        
        # 为每个阳离子和目标阴离子生成RDF对
        for idx, element in enumerate(element_list, start=1):
            if element in cation_elements:
                for anion_idx, anion_element in enumerate(element_list, start=1):
                    if anion_element in ['O', 'N', 'P', 'B']:
                        rdf_pairs += f'{idx} {anion_idx} '
                        rdf_labels += f'rdf_{anion_element} '
        
        return element_list, rdf_pairs.strip(), rdf_labels.strip()
    
    def create_lammps_in_list(self, destination_path, system_name, element_list, rdf_pairs):
        """创建LAMMPS .in.list文件
        
        Args:
            destination_path: 目标目录路径
            system_name: 系统名称
            element_list: 元素列表
            rdf_pairs: RDF对字符串
            
        Returns:
            .in.list文件路径
        """
        in_list_content = f'variable element_list index "{" ".join(element_list)}"\n'
        in_list_content += f'variable rdf_pair string "{rdf_pairs}"\n'
        
        in_list_path = destination_path / f"{system_name}.in.list"
        with open(in_list_path, 'w', encoding='utf-8') as in_list_file:
            in_list_file.write(in_list_content)
            
        self.logger.info(f"Generated LAMMPS .in.list file: {in_list_path}")
        return in_list_path
    
    def create_lammps_input(self, system, destination_path, element_list, rdf_pairs, rdf_labels, simulation_params):
        """创建LAMMPS输入文件
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            element_list: 元素列表
            rdf_pairs: RDF对字符串
            rdf_labels: RDF标签字符串
            simulation_params: 模拟参数字典
            
        Returns:
            LAMMPS输入文件路径
        """
        lmp_content = "# ----------------- Variable Section -----------------\n"
        
        # 定义分子变量
        for i, cation in enumerate(system.cations, start=1):
            if cation.is_valid():
                lmp_content += f'variable Cation{i} index {cation.name}\n'
                
        for i, anion in enumerate(system.anions, start=1):
            if anion.is_valid():
                lmp_content += f'variable Anion{i} index {anion.name}\n'
                
        for i, solvent in enumerate(system.solvents, start=1):
            if solvent.is_valid():
                lmp_content += f'variable Solvent{i} index {solvent.name}\n'
        
        # 定义文件变量
        lmp_content += f"""
variable infile string {system.name}
variable outname string {system.name}

include {system.name}.in.init
read_data {system.name}.data
include {system.name}.in.settings
include {system.name}.in.list
include {system.name}.in.list_salt

variable Nsteps_NPT equal {simulation_params["npt"]["steps"]}
variable Nsteps_NVT equal {simulation_params["nvt"]["steps"]}

variable Freq_trj_npt index {simulation_params["npt"]["freq_trj"]}
variable Freq_trj_nvt index {simulation_params["nvt"]["freq_trj"]}

variable Temp_NPT equal {system.temperature}
variable Temp_NVT equal {system.temperature}

variable thermo_freq equal {simulation_params["general"]["thermo_freq"]}
variable timestep equal {simulation_params["general"]["timestep"]}
"""
        
        # 定义分子组
        for i, cation in enumerate(system.cations, start=1):
            if cation.is_valid():
                lmp_content += f'group Cation{i} union {cation.name}\n'
                
        for i, anion in enumerate(system.anions, start=1):
            if anion.is_valid():
                lmp_content += f'group Anion{i} union {anion.name}\n'
        
        # 热力学设置
        lmp_content += f"""
thermo_style custom step cpu cpuremain temp density lx ly lz etotal ke pe evdwl ecoul elong ebond eangle edihed eimp
thermo ${{thermo_freq}}
timestep ${{timestep}}

minimize {simulation_params["minimization"]["etol"]} {simulation_params["minimization"]["ftol"]} {simulation_params["minimization"]["maxiter"]} {simulation_params["minimization"]["maxeval"]}

write_data ${{infile}}_after_minimize.data nocoeff
write_dump all custom ${{infile}}_after_minimize.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id

reset_timestep 0

dump trj_npt all custom ${{Freq_trj_npt}} NPT_${{outname}}.lammpstrj id element mol type x y z q
dump_modify trj_npt flush yes element ${{element_list}} sort id

fix fxnpt all npt temp ${{Temp_NPT}} ${{Temp_NPT}} $(100.0*dt) iso 0.0 0.0 $(1000.0*dt)
run ${{Nsteps_NPT}}
unfix fxnpt

undump trj_npt

write_data ${{infile}}_after_npt.data
write_restart ${{infile}}_restart_after_npt.data
write_dump all custom ${{infile}}_after_npt.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id

reset_timestep 0
"""
        
        # MSD计算设置
        for i, anion in enumerate(system.anions, start=1):
            if anion.is_valid():
                lmp_content += f'compute An{i} {anion.name} msd com yes\n'
                lmp_content += f'fix An{i}msd {anion.name} ave/time {simulation_params["nvt"]["freq_trj"]} 1 {simulation_params["nvt"]["freq_trj"]} c_An{i}[1] c_An{i}[2] c_An{i}[3] c_An{i}[4] file out_{anion.name}_msd.dat title1 "t msd msd msd msd_{anion.name}" title2 "fs {anion.name}_x {anion.name}_y {anion.name}_z {anion.name}_total"\n'
                
        for i, cation in enumerate(system.cations, start=1):
            if cation.is_valid():
                lmp_content += f'compute Ca{i} {cation.name} msd com yes\n'
                lmp_content += f'fix Ca{i}msd {cation.name} ave/time {simulation_params["nvt"]["freq_trj"]} 1 {simulation_params["nvt"]["freq_trj"]} c_Ca{i}[1] c_Ca{i}[2] c_Ca{i}[3] c_Ca{i}[4] file out_{cation.name}_msd.dat title1 "t msd msd msd msd_{cation.name}" title2 "fs {cation.name}_x {cation.name}_y {cation.name}_z {cation.name}_total"\n'
        
        # RDF计算和NVT模拟设置
        lmp_content += f"""
compute rdfc1 all rdf 100 ${{rdf_pair}}
fix rdff1 all ave/time $(v_Nsteps_NVT/1000) 1000 ${{Nsteps_NVT}} c_rdfc1[*] file out_rdf.dat mode vector title3 "RDF {rdf_labels}"

dump trj_nvt all custom ${{Freq_trj_nvt}} NVT_${{outname}}.lammpstrj id element mol type x y z q
dump_modify trj_nvt flush yes element ${{element_list}} sort id
dump utrj_nvt all custom ${{Freq_trj_nvt}} NVT_${{outname}}_un.lammpstrj id element mol type xu yu zu ix iy iz q
dump_modify utrj_nvt flush yes element ${{element_list}} sort id

fix fxnvt all nvt temp ${{Temp_NVT}} ${{Temp_NVT}} $(100.0*dt)
run ${{Nsteps_NVT}}
unfix fxnvt

undump trj_nvt
undump utrj_nvt

write_data ${{infile}}_after_nvt.data
write_restart ${{infile}}_restart_after_nvt.data
write_dump all custom ${{infile}}_after_nvt.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id
"""
        
        lmp_input_path = destination_path / f"{system.name}.in"
        with open(lmp_input_path, 'w') as lmp_file:
            lmp_file.write(lmp_content)
            
        self.logger.info(f"Generated LAMMPS input file: {lmp_input_path}")
        return lmp_input_path

class LAMMPSFileGenerator:
    """LAMMPS文件生成器类，专门处理LAMMPS相关的输入文件生成"""
    
    def __init__(self):
        """初始化LAMMPS文件生成器"""
        self.logger = Logger().get_logger()
    
    def generate_input_file(self, config, output_path):
        """生成LAMMPS输入文件
        
        Args:
            config: 配置数据，包含电解质配方等信息
            output_path: 输出文件路径
            
        Returns:
            生成的输入文件路径
        """
        self.logger.info(f"生成LAMMPS输入文件: {output_path}")
        
        formulation_name = config.get('formulation_name', 'Electrolyte_Formulation')
        salts = config.get('salts', [])
        solvents = config.get('solvents', [])
        temperature = config.get('temperature', 298)
        pressure = config.get('pressure', 1.0)
        time_step = config.get('time_step', 2.0)
        equilibration_steps = config.get('equilibration_steps', 500000)
        production_steps = config.get('production_steps', 1000000)
        cutoff = config.get('cutoff', 12.0)
        
        # 创建输出目录
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # 生成输入文件内容
        content = [
            f"# LAMMPS输入文件 - {formulation_name}",
            f"# 由Molyte-Cursor生成\n",
            "units real",
            "atom_style full",
            f"timestep {time_step}",
            f"pair_style lj/cut/coul/long {cutoff}\n",
            "# 溶剂组分:"
        ]
        
        for solvent in solvents:
            content.append(f"# - {solvent['name']}: {solvent['concentration']}")
        
        content.append("\n# 盐组分:")
        for salt in salts:
            content.append(f"# - {salt['name']}: {salt['concentration']} mol/L")
        
        content.append(f"\n# 模拟参数")
        content.append(f"# 温度: {temperature} K")
        content.append(f"# 压力: {pressure} atm")
        content.append(f"# 平衡步数: {equilibration_steps}")
        content.append(f"# 生产步数: {production_steps}")
        content.append(f"# 截断距离: {cutoff} Å")
        
        # 写入文件
        with open(output_path, 'w') as f:
            f.write('\n'.join(content))
        
        self.logger.info(f"成功生成LAMMPS输入文件: {output_path}")
        return output_path
        
    def generate_inp_file(self, config, output_path):
        """生成INP格式的电解质输入文件
        
        Args:
            config: 配置数据，包含电解质配方等信息
            output_path: 输出文件路径
            
        Returns:
            生成的INP文件路径
        """
        self.logger.info(f"生成INP输入文件: {output_path}")
        
        # 提取配置信息
        formulation_name = config.get('formulation_name', 'Electrolyte_Formulation')
        salts = config.get('salts', [])
        solvents = config.get('solvents', [])
        temperature = config.get('temperature', 298)
        pressure = config.get('pressure', 1.0)
        time_step = config.get('time_step', 2.0)
        equilibration_steps = config.get('equilibration_steps', 500000)
        production_steps = config.get('production_steps', 1000000)
        cutoff = config.get('cutoff', 12.0)
        box_size = config.get('box_size', 50)
        concentration = config.get('concentration', 1.0)
        
        # 创建输出目录
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # 生成INP文件内容
        content = [
            "START",
            f"formulaName: {formulation_name}",
            f"temperature: {temperature}",
            f"Box_size: {box_size}",
            f"concentration: {concentration}",
            f"time_step: {time_step}",
            f"equilibration_steps: {equilibration_steps}",
            f"production_steps: {production_steps}",
            f"cutoff: {cutoff}",
            f"pressure: {pressure}",
            "",
            "cations:"
        ]
        
        # 添加阳离子
        if salts:
            for salt in salts:
                content.append(f"  - name: {salt.get('cation', 'Li')}")
                content.append(f"    charge: 1")
                content.append(f"    concentration: {salt.get('concentration', 1.0)}")
        else:
            content.append("  - name: Li")
            content.append("    charge: 1")
            content.append("    concentration: 1.0")
        
        # 添加阴离子
        content.append("")
        content.append("anions:")
        if salts:
            for salt in salts:
                content.append(f"  - name: {salt.get('anion', 'PF6')}")
                content.append(f"    charge: -1")
                content.append(f"    concentration: {salt.get('concentration', 1.0)}")
        else:
            content.append("  - name: PF6")
            content.append("    charge: -1")
            content.append("    concentration: 1.0")
        
        # 添加溶剂
        content.append("")
        content.append("solvents:")
        if solvents:
            for solvent in solvents:
                content.append(f"  - name: {solvent.get('name', 'EC')}")
                content.append(f"    smile: {solvent.get('smile', 'C1OC(=O)O1')}")
                content.append(f"    ratio: {solvent.get('concentration', 1.0)}")
        else:
            content.append("  - name: EC")
            content.append("    smile: C1OC(=O)O1")
            content.append("    ratio: 1.0")
        
        content.append("")
        content.append("END")
        
        # 写入文件
        with open(output_path, 'w') as f:
            f.write('\n'.join(content))
        
        self.logger.info(f"成功生成INP输入文件: {output_path}")
        return output_path

# 添加模块级别的函数，用于从Django后端调用
def generate_input_files(config, output_dir):
    """生成输入文件的函数
    
    Args:
        config: 配置数据，包含电解质配方等信息
        output_dir: 输出目录路径
        
    Returns:
        生成的文件路径字典
    """
    logger = Logger().get_logger()
    logger.info(f"开始生成输入文件，输出目录: {output_dir}")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 初始化文件生成器
    lammps_generator = LAMMPSFileGenerator()
    
    # 生成INP文件
    inp_path = os.path.join(output_dir, "input.inp")
    inp_file = lammps_generator.generate_inp_file(config, inp_path)
    
    # 生成LAMMPS输入文件
    lammps_path = os.path.join(output_dir, "input.lammps")
    lammps_file = lammps_generator.generate_input_file(config, lammps_path)
    
    # 返回生成的文件路径
    return {
        "inp_file": inp_file,
        "lammps_file": lammps_file
    } 