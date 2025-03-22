import os
import logging


class LAMMPSGenerator:
    """LAMMPS输入文件生成器类"""
    
    def __init__(self):
        """初始化LAMMPS生成器"""
        self.logger = logging.getLogger("LAMMPSGenerator")
    
    def generate_lammps_input(self, system_name, box_size, output_dir):
        """生成LAMMPS输入文件
        
        Args:
            system_name: 系统名称
            box_size: 模拟盒子大小
            output_dir: 输出目录
            
        Returns:
            生成的LAMMPS输入文件路径
        """
        self.logger.info(f"为系统 {system_name} 生成LAMMPS输入文件")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 准备LAMMPS输入文件路径
        lammps_file = os.path.join(output_dir, f"{system_name}.in")
        
        # 创建LAMMPS输入文件
        with open(lammps_file, 'w') as f:
            f.write(f"# LAMMPS input file for {system_name}\n\n")
            
            # 基本设置
            f.write("# 初始化设置\n")
            f.write("units real\n")
            f.write("atom_style full\n")
            f.write("bond_style harmonic\n")
            f.write("angle_style harmonic\n")
            f.write("dihedral_style opls\n")
            f.write("improper_style harmonic\n")
            f.write("pair_style lj/cut/coul/long 10.0\n")
            f.write("kspace_style pppm 1.0e-4\n\n")
            
            # 读取数据文件
            f.write("# 读取数据文件\n")
            f.write(f"read_data {system_name}.data\n\n")
            
            # 分子设置
            f.write("# 设置分子\n")
            f.write("group solute type 1 2\n")
            f.write("group solvent type 3 4 5\n\n")
            
            # 最小化
            f.write("# 能量最小化\n")
            f.write("minimize 1.0e-4 1.0e-6 1000 10000\n\n")
            
            # 平衡设置
            f.write("# 设置平衡运行\n")
            f.write("timestep 1.0\n")
            f.write("fix 1 all nvt temp 300.0 300.0 100.0\n")
            f.write("thermo 1000\n")
            f.write("thermo_style custom step temp press pe ke etotal density\n\n")
            
            # 运行平衡
            f.write("# 运行平衡\n")
            f.write("run 10000\n\n")
            
            # 生产运行
            f.write("# 生产运行\n")
            f.write("reset_timestep 0\n")
            f.write("dump 1 all custom 1000 traj.lammpstrj id mol type x y z vx vy vz\n")
            f.write("run 50000\n")
        
        self.logger.info(f"已生成LAMMPS输入文件: {lammps_file}")
        return lammps_file 