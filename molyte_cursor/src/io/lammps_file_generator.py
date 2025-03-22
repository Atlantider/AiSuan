import os
import logging
from .lammps_generator import LAMMPSGenerator

class LAMMPSFileGenerator(LAMMPSGenerator):
    """LAMMPS输入文件生成器类，扩展自LAMMPSGenerator"""
    
    def __init__(self):
        """初始化LAMMPS文件生成器"""
        super().__init__()
        self.logger = logging.getLogger("LAMMPSFileGenerator")
    
    def generate_input_file(self, config, output_path, data_file=None, topology_files=None, forcefield_files=None):
        """生成LAMMPS输入文件
        
        Args:
            config: 配置数据，包含电解质配方等信息
            output_path: 输出文件路径
            data_file: Packmol生成的LAMMPS数据文件路径
            topology_files: 分子拓扑文件列表
            forcefield_files: 力场参数文件列表
            
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
        box_size = config.get('Box_size', 50.0)
        
        # 创建输出目录
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # 生成输入文件内容
        content = [
            f"# LAMMPS输入文件 - {formulation_name}",
            f"# 由Molyte-Cursor生成\n",
            "# 基本设置",
            "units real",
            "atom_style full",
            "dimension 3",
            "boundary p p p",
            "dielectric 1.0",
            f"timestep {time_step}",
            f"pair_style lj/cut/coul/long {cutoff}",
            "special_bonds lj/coul 0.0 0.0 0.5"
        ]
        
        # 添加数据文件读取
        if data_file and os.path.exists(data_file):
            data_file_rel = os.path.relpath(data_file, os.path.dirname(output_path))
            content.append(f"\n# 读取系统数据文件")
            content.append(f"read_data {data_file_rel}")
        else:
            content.append(f"\n# 警告：未提供数据文件")
            
        # 添加力场参数
        if forcefield_files:
            content.append(f"\n# 读取力场参数")
            for ff_file in forcefield_files:
                if os.path.exists(ff_file):
                    ff_file_rel = os.path.relpath(ff_file, os.path.dirname(output_path))
                    content.append(f"include {ff_file_rel}")
        
        # 添加分子组定义 - 根据配置动态确定
        content.append("\n# 定义分子组")
        
        # 从配置中获取分子类型信息
        cation_types = []
        anion_types = []
        solvent_types = []
        all_types = set()
        
        # 尝试从配置中获取阳离子类型
        for i in range(1, 10):  # 最多10个阳离子
            cation_key = f"Cation{i}_name"
            if cation_key in config and config[cation_key]:
                # 阳离子的类型基于索引，从1开始
                cation_type = i
                cation_types.append(cation_type)
                all_types.add(cation_type)
                
        # 尝试从配置中获取阴离子类型
        for i in range(1, 10):  # 最多10个阴离子
            anion_key = f"Anion{i}_name"
            if anion_key in config and config[anion_key]:
                # 阴离子的类型基于索引，跟在阳离子后面
                anion_type = len(cation_types) + i
                anion_types.append(anion_type)
                all_types.add(anion_type)
                
        # 尝试从配置中获取溶剂类型
        for i in range(1, 10):  # 最多10个溶剂
            solvent_key = f"Sol{i}_name"
            if solvent_key in config and config[solvent_key]:
                # 溶剂的类型基于索引，跟在阴离子后面
                solvent_type = len(cation_types) + len(anion_types) + i
                solvent_types.append(solvent_type)
                all_types.add(solvent_type)
        
        # 如果配置中没有分子类型信息，使用默认值
        if not all_types:
            content.append("group cations type 1 2")
            content.append("group anions type 3 4")
            content.append("group ions union cations anions")
            content.append("group solvents type 5 6 7")
        else:
            if cation_types:
                content.append(f"group cations type {' '.join(map(str, cation_types))}")
            if anion_types:
                content.append(f"group anions type {' '.join(map(str, anion_types))}")
            if cation_types or anion_types:
                content.append("group ions union cations anions")
            if solvent_types:
                content.append(f"group solvents type {' '.join(map(str, solvent_types))}")
        
        # 添加配对样式和力场参数
        content.append("\n# 设置相互作用")
        content.append("kspace_style pppm 1.0e-4")
        
        # 最小化能量
        content.append("\n# 能量最小化")
        content.append("minimize 1.0e-4 1.0e-6 1000 10000")
        
        # 设置温度和压力控制
        content.append("\n# 设置温度和压力控制")
        content.append(f"variable T equal {temperature}")
        content.append("variable dt equal 2.0")
        content.append("variable tdamp equal 100*${dt}")
        content.append("variable pdamp equal 1000*${dt}")
        
        # 平衡阶段
        content.append("\n# 平衡阶段 - NVT")
        content.append("fix 1 all nvt temp ${T} ${T} ${tdamp}")
        content.append("thermo 1000")
        content.append("thermo_style custom step temp press density etotal")
        content.append(f"run {equilibration_steps // 2}")
        content.append("unfix 1")
        
        content.append("\n# 平衡阶段 - NPT")
        content.append(f"fix 1 all npt temp ${{T}} ${{T}} ${{tdamp}} iso {pressure} {pressure} ${{pdamp}}")
        content.append(f"run {equilibration_steps // 2}")
        
        # 生产阶段
        content.append("\n# 生产阶段 - NPT")
        content.append("reset_timestep 0")
        
        # 计算MSD
        content.append("\n# 计算扩散系数")
        content.append("compute msd_cations cations msd")
        content.append("compute msd_anions anions msd")
        content.append("compute msd_solvents solvents msd")
        content.append("compute msd_ions ions msd")
        content.append("compute msd_all all msd")
        
        # 记录轨迹和热力学信息
        content.append("\n# 记录轨迹和热力学信息")
        content.append("dump 1 all custom 5000 traj.lammpstrj id mol type x y z vx vy vz")
        content.append("dump_modify 1 sort id")
        content.append("thermo 5000")
        content.append("thermo_style custom step temp press density etotal c_msd_cations[4] c_msd_anions[4] c_msd_solvents[4] c_msd_all[4]")
        
        # 记录MSD数据
        content.append("fix msd_output all print 5000 \"${time} ${temp} ${press} ${density} ${etotal} ${c_msd_cations[4]} ${c_msd_anions[4]} ${c_msd_solvents[4]} ${c_msd_all[4]}\" file msd.dat screen no")
        
        # 运行生产阶段
        content.append(f"\nrun {production_steps}")
        
        # 写入输入文件
        with open(output_path, 'w') as f:
            f.write('\n'.join(content))
        
        self.logger.info(f"成功生成LAMMPS输入文件: {output_path}")
        return output_path

    def generate_data_file(self, system_info, output_path):
        """生成LAMMPS数据文件
        
        Args:
            system_info: 系统信息字典
            output_path: 输出文件路径
            
        Returns:
            生成的数据文件路径
        """
        self.logger.info(f"生成LAMMPS数据文件: {output_path}")
        
        # 这是一个简单的实现，实际应用中可能需要更复杂的逻辑
        with open(output_path, 'w') as f:
            f.write(f"# LAMMPS data file for {system_info.get('name', 'system')}\n\n")
            f.write(f"{system_info.get('atoms', 0)} atoms\n")
            f.write(f"{system_info.get('atom_types', 0)} atom types\n\n")
            
            # 添加盒子尺寸信息
            box = system_info.get('box', {'xlo': 0, 'xhi': 10, 'ylo': 0, 'yhi': 10, 'zlo': 0, 'zhi': 10})
            f.write(f"{box.get('xlo', 0)} {box.get('xhi', 10)} xlo xhi\n")
            f.write(f"{box.get('ylo', 0)} {box.get('yhi', 10)} ylo yhi\n")
            f.write(f"{box.get('zlo', 0)} {box.get('zhi', 10)} zlo zhi\n\n")
            
            # 添加更多数据内容...
            
        return output_path
    
    def combine_files(self, input_files, data_file, output_file):
        """将多个输入文件合并为一个LAMMPS脚本
        
        Args:
            input_files: 输入文件列表
            data_file: 数据文件路径
            output_file: 输出文件路径
            
        Returns:
            合并后的文件路径
        """
        self.logger.info(f"合并LAMMPS文件: {output_file}")
        
        with open(output_file, 'w') as out:
            out.write(f"# Combined LAMMPS input script\n\n")
            out.write(f"read_data {os.path.basename(data_file)}\n\n")
            
            # 添加其他输入文件的内容
            for in_file in input_files:
                if os.path.exists(in_file):
                    out.write(f"# Including file: {os.path.basename(in_file)}\n")
                    with open(in_file, 'r') as f:
                        out.write(f.read())
                    out.write("\n\n")
            
        return output_file 