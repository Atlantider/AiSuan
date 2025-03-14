"""
模拟处理器模块，负责执行分子模拟的全过程
"""
import os
import subprocess
import json
import hashlib
from pathlib import Path
from ..utils.logger import Logger
from ..utils.command_executor import CommandExecutor
from ..utils.config_loader import ConfigLoader
from ..utils.charge_modifier import ChargeModifier
from ..io.file_generator import FileGenerator

class Simulator:
    """模拟处理器类"""
    
    # 支持的阳离子和阴离子列表
    SUPPORTED_CATIONS = ['Li', 'Na', 'K', 'Zn', 'Mg', 'Mg0.8', 'Zn0.8']
    SUPPORTED_ANIONS = ['PF6', 'FSI', 'TFSI', 'DFOB', 'Cl', 'NO3', 'ClO4', 'BF4']
    
    def __init__(self):
        """初始化模拟处理器"""
        self.logger = Logger().get_logger()
        self.executor = CommandExecutor()
        self.config = ConfigLoader()
        self.charge_modifier = ChargeModifier()
        self.file_generator = FileGenerator()
        
        # 加载配置
        self.config.load_paths_config()
        self.config.load_simulation_params()
        
        # 获取路径配置
        self.file_generate_path = self.config.get_path('file_generate')
        self.initial_salts_path = self.config.get_path('initial_salts')
        self.ligpargen_path = self.config.get_path('ligpargen')
        self.packmol_path = self.config.get_path('packmol')
        self.lempyfy_path = self.config.get_path('lempyfy')
        self.moltemplate_path = self.config.get_path('moltemplate')
        self.resp_path = self.config.get_path('resp')
        self.charge_save_path = self.config.get_path('charge_save')
        
        # 创建必要的目录
        self.file_generate_path.mkdir(parents=True, exist_ok=True)
        self.charge_save_path.mkdir(parents=True, exist_ok=True)
        
        # 加载计算历史
        self.calculation_history_file = Path("calculation_history.json")
        self.calculation_history = self._load_calculation_history()
        
        # 初始化拒绝组分列表
        self.rejected_components = []
    
    def _load_calculation_history(self):
        """加载计算历史记录
        
        Returns:
            计算历史字典
        """
        if self.calculation_history_file.exists():
            with open(self.calculation_history_file, 'r') as f:
                try:
                    return json.load(f)
                except json.JSONDecodeError:
                    self.logger.error(f"无法解析计算历史文件: {self.calculation_history_file}")
                    return {}
        return {}
    
    def save_calculation_history(self, system_hash, path):
        """保存计算历史记录
        
        Args:
            system_hash: 系统哈希值
            path: 计算路径
        """
        self.calculation_history[system_hash] = str(path)
        with open(self.calculation_history_file, 'w') as f:
            json.dump(self.calculation_history, f, indent=2)
            self.logger.info(f"计算历史已保存到: {self.calculation_history_file}")
    
    def save_rejected_components(self):
        """保存被拒绝的组分列表"""
        with open("rejected_components.txt", "a") as f:
            for item in self.rejected_components:
                f.write(f"{item['system_name']} - {item['component_type']}: {item['component_name']} - {item['reason']}\n")
        self.logger.info("已更新拒绝组分列表")
    
    def is_supported_ion(self, name, ion_type):
        """检查离子是否被支持
        
        Args:
            name: 离子名称
            ion_type: 离子类型('cation'或'anion')
            
        Returns:
            布尔值，True表示支持
        """
        if ion_type == 'cation':
            return name in self.SUPPORTED_CATIONS
        elif ion_type == 'anion':
            return name in self.SUPPORTED_ANIONS
        return True  # 溶剂默认支持
    
    def prepare_destination(self, system):
        """准备目标目录
        
        Args:
            system: 分子系统对象
            
        Returns:
            目标目录路径
        """
        # 在目录名中添加温度信息
        temperature = getattr(system, 'temperature', 300.0)
        folder_name = f"{system.name}_{int(temperature)}K"
        destination_path = self.file_generate_path / folder_name
        
        # 如果目录已存在，则删除
        if destination_path.exists():
            self.logger.info(f"Directory {destination_path} already exists. Removing it.")
            self.executor.run_command(f'rm -rf {destination_path}')
        
        # 复制初始盐文件
        self.executor.run_command(f'cp -rf {self.initial_salts_path} {self.file_generate_path}')
        self.executor.run_command(f'mv -f {self.file_generate_path}/inital_salts/ {destination_path}')
        
        self.logger.info(f"Prepared destination directory: {destination_path}")
        return destination_path
    
    def validate_components(self, system):
        """验证系统组分是否都被支持
        
        Args:
            system: 分子系统对象
            
        Returns:
            (布尔值, 拒绝列表) - True表示所有组分都被支持
        """
        all_supported = True
        rejected = []
        
        # 检查阳离子
        for cation in system.cations:
            if not self.is_supported_ion(cation.name, 'cation'):
                all_supported = False
                rejected.append({
                    'system_name': system.name,
                    'component_type': 'cation',
                    'component_name': cation.name,
                    'reason': '不支持的阳离子类型'
                })
                self.logger.warning(f"不支持的阳离子: {cation.name}，将被记录但继续处理")
        
        # 检查阴离子
        for anion in system.anions:
            if not self.is_supported_ion(anion.name, 'anion'):
                all_supported = False
                rejected.append({
                    'system_name': system.name,
                    'component_type': 'anion',
                    'component_name': anion.name,
                    'reason': '不支持的阴离子类型'
                })
                self.logger.warning(f"不支持的阴离子: {anion.name}，将被记录但继续处理")
        
        # 检查盒子尺寸
        box_size = getattr(system, 'box_size', 40.0)
        if isinstance(box_size, (int, float)) and (box_size <= 0 or box_size > 100):
            all_supported = False
            rejected.append({
                'system_name': system.name,
                'component_type': 'system_parameter',
                'component_name': 'box_size',
                'reason': f"无效的盒子尺寸: {box_size}，必须在0-100范围内"
            })
            self.logger.warning(f"盒子尺寸 {box_size} 不在有效范围内(0-100埃)，将被记录但继续处理")
        
        return all_supported, rejected
    
    def prepare_solvents(self, system, destination_path):
        """准备溶剂分子
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
        """
        # 处理每个溶剂
        for solvent in system.solvents:
            if solvent.is_valid() and hasattr(solvent, 'smile') and solvent.smile:
                charge_file_path = self.charge_save_path / f"{solvent.name}.charmm.chg"
                
                # 使用ligpargen生成初始电荷文件
                self.executor.run_command(
                    f"{self.ligpargen_path}/ligpargen -s '{solvent.smile}' -n {solvent.name} -r MOL -c 0 -o 0 -cgen CM1A",
                    cwd=str(destination_path)
                )
                
                # 如果电荷文件不存在，生成电荷文件
                if not charge_file_path.exists():
                    self.logger.info(f"Generating charges for {solvent.name}")
                    self.executor.run_command(
                        f"{self.resp_path}/RESP2.sh '{solvent.name}.charmm.pdb'",
                        cwd=str(destination_path)
                    )
                    self.executor.run_command(
                        f"cp {solvent.name}.charmm.chg {charge_file_path}",
                        cwd=str(destination_path)
                    )
                else:
                    self.logger.info(f"Charge file for {solvent.name} already exists. Skipping generation.")
                
                # 修改LAMMPS文件中的电荷
                try:
                    self.charge_modifier.modify_lmp_charges(
                        str(destination_path / f"{solvent.name}.lammps.lmp"),
                        charge_file_path
                    )
                except Exception as e:
                    self.logger.error(f"Failed to modify charges for {solvent.name}: {str(e)}")
            else:
                self.logger.warning(f"溶剂 {solvent.name} 缺少SMILE字符串，跳过电荷生成")
                self.rejected_components.append({
                    'system_name': system.name,
                    'component_type': 'solvent',
                    'component_name': solvent.name,
                    'reason': '缺少SMILE字符串'
                })
    
    def run_packmol(self, system, destination_path):
        """运行Packmol程序
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            
        Returns:
            是否成功
        """
        # 创建Packmol输入文件
        packmol_input_file = self.file_generator.create_packmol_input(system, destination_path)
        
        # 运行Packmol
        try:
            self.executor.run_command(f'{self.packmol_path} < {packmol_input_file}', cwd=str(destination_path))
            self.logger.info("Packmol completed successfully")
            return True
        except Exception as e:
            self.logger.error(f"Packmol failed: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'packmol',
                'reason': f"Packmol运行失败: {str(e)}"
            })
            return False
    
    def prepare_moltemplate_files(self, system, destination_path):
        """准备Moltemplate文件
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            
        Returns:
            是否成功
        """
        # 为每个溶剂生成Moltemplate文件
        for solvent in system.solvents:
            if solvent.is_valid():
                try:
                    self.executor.run_command(
                        f"python {self.lempyfy_path} {solvent.name}.lammps.lmp > {solvent.name}.lt",
                        cwd=str(destination_path)
                    )
                except Exception as e:
                    self.logger.error(f"Failed to generate LT file for {solvent.name}: {str(e)}")
                    self.rejected_components.append({
                        'system_name': system.name,
                        'component_type': 'process',
                        'component_name': f'lt_generation_{solvent.name}',
                        'reason': f"生成LT文件失败: {str(e)}"
                    })
        
        # 创建主Moltemplate文件
        main_lt_file = self.file_generator.create_main_lt_file(system, destination_path)
        
        # 运行Moltemplate
        try:
            self.executor.run_command(
                f"{self.moltemplate_path} {main_lt_file}",
                cwd=str(destination_path)
            )
            self.logger.info("Moltemplate completed successfully")
            return True
        except Exception as e:
            self.logger.error(f"Moltemplate failed: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'moltemplate',
                'reason': f"Moltemplate运行失败: {str(e)}"
            })
            return False
    
    def prepare_lammps_files(self, system, destination_path):
        """准备LAMMPS文件
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            
        Returns:
            是否成功
        """
        # 检查系统LAMMPS数据文件是否存在
        data_file = destination_path / "system.data"
        if not data_file.exists():
            self.logger.error(f"LAMMPS data file not found: {data_file}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'file',
                'component_name': 'system.data',
                'reason': "LAMMPS数据文件不存在"
            })
            return False
        
        # 提取元素列表
        try:
            with open(destination_path / "log.cite", "r") as f:
                element_list = []
                for line in f:
                    if "Atoms in final file:" in line:
                        element_str = line.split(":")[1].strip()
                        element_list = [e.strip() for e in element_str.split(',')]
                        break
            self.logger.info(f"Element list: {element_list}")
        except Exception as e:
            self.logger.error(f"Failed to extract element list: {str(e)}")
            element_list = []
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'element_extraction',
                'reason': f"提取元素列表失败: {str(e)}"
            })
        
        # 解析RDF配对
        rdf_pairs = []
        rdf_labels = []
        try:
            in_list_file = destination_path / "Atoms_In_List.txt"
            if in_list_file.exists():
                with open(in_list_file, 'r') as f:
                    for line in f:
                        if line.startswith('RDF_PAIR'):
                            parts = line.strip().split(':')
                            if len(parts) >= 2:
                                pair_parts = parts[1].strip().split('-')
                                if len(pair_parts) == 2:
                                    rdf_pairs.append((pair_parts[0].strip(), pair_parts[1].strip()))
                                    rdf_labels.append(parts[0].strip())
            self.logger.info(f"RDF pairs: {rdf_pairs}")
            self.logger.info(f"RDF labels: {rdf_labels}")
        except Exception as e:
            self.logger.error(f"Failed to parse RDF pairs: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'rdf_extraction',
                'reason': f"解析RDF配对失败: {str(e)}"
            })
        
        # 准备模拟参数
        temperature = getattr(system, 'temperature', 300.0)
        simulation_params = {
            "temperature": temperature,
            "equilibration": {
                "timestep": self.config.get_simulation_param("equilibration.timestep"),
                "duration": self.config.get_simulation_param("equilibration.duration"),
                "dump_freq": self.config.get_simulation_param("equilibration.dump_freq")
            },
            "production": {
                "timestep": self.config.get_simulation_param("production.timestep"),
                "duration": self.config.get_simulation_param("production.duration"),
                "dump_freq": self.config.get_simulation_param("production.dump_freq"),
                "thermo_freq": self.config.get_simulation_param("production.thermo_freq")
            },
            "minimization": {
                "etol": self.config.get_simulation_param("minimization.etol"),
                "ftol": self.config.get_simulation_param("minimization.ftol"),
                "maxiter": self.config.get_simulation_param("minimization.maxiter"),
                "maxeval": self.config.get_simulation_param("minimization.maxeval")
            }
        }
        
        # 创建LAMMPS输入文件
        try:
            self.file_generator.create_lammps_input(
                system, destination_path, element_list, rdf_pairs, rdf_labels, simulation_params
            )
            return True
        except Exception as e:
            self.logger.error(f"Failed to create LAMMPS input files: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'lammps_input',
                'reason': f"创建LAMMPS输入文件失败: {str(e)}"
            })
            return False
    
    def submit_job(self, system, destination_path, submit=False):
        """提交计算作业
        
        Args:
            system: 分子系统对象
            destination_path: 目标目录路径
            submit: 是否实际提交作业
            
        Returns:
            是否成功
        """
        # 创建作业脚本
        try:
            job_script_path = self.file_generator.create_job_script(destination_path, system.name)
        except Exception as e:
            self.logger.error(f"Failed to create job script: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'job_script',
                'reason': f"创建作业脚本失败: {str(e)}"
            })
            return False
        
        # 如果需要提交作业
        if submit:
            try:
                self.executor.run_command(f'qsub {job_script_path}', cwd=str(destination_path))
                self.logger.info(f"Submitted job for {system.name}")
                return True
            except Exception as e:
                self.logger.error(f"Failed to submit job: {str(e)}")
                self.rejected_components.append({
                    'system_name': system.name,
                    'component_type': 'process',
                    'component_name': 'job_submission',
                    'reason': f"提交作业失败: {str(e)}"
                })
                return False
        else:
            self.logger.info(f"Job script prepared but not submitted for {system.name}")
            return True
    
    def generate_system_hash(self, system):
        """生成系统唯一哈希值
        
        Args:
            system: 分子系统对象
            
        Returns:
            哈希字符串
        """
        # 提取用于哈希的关键数据
        hash_data = {
            'temperature': getattr(system, 'temperature', 300.0),
            'components': []
        }
        
        # 添加所有组分信息
        for component_list, component_type in [
            (system.cations, 'cation'),
            (system.anions, 'anion'),
            (system.solvents, 'solvent')
        ]:
            for comp in component_list:
                if comp.is_valid():
                    hash_data['components'].append({
                        'name': comp.name,
                        'smile': getattr(comp, 'smile', ''),
                        'ratio': getattr(comp, 'ratio', 1.0),
                        'type': component_type
                    })
        
        # 生成哈希值
        hash_str = json.dumps(hash_data, sort_keys=True)
        return hashlib.md5(hash_str.encode()).hexdigest()
    
    def check_duplicate_calculation(self, system):
        """检查是否已经计算过相同的系统
        
        Args:
            system: 分子系统对象
            
        Returns:
            (是否是重复计算, 先前计算的路径)
        """
        system_hash = self.generate_system_hash(system)
        
        if system_hash in self.calculation_history:
            previous_path = self.calculation_history[system_hash]
            self.logger.warning(f"系统已经计算过，路径: {previous_path}")
            return True, previous_path, system_hash
        
        return False, None, system_hash
    
    def simulate(self, system, submit_job=False):
        """执行完整的模拟流程
        
        Args:
            system: 分子系统对象
            submit_job: 是否提交作业
            
        Returns:
            是否成功
        """
        self.logger.info(f"Starting simulation for system: {system.name}")
        
        # 检查系统是否有效
        if not system.has_valid_cation():
            self.logger.warning(f"System {system.name} has no valid cation. Skipping.")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'system',
                'component_name': system.name,
                'reason': "没有有效的阳离子"
            })
            return False
        
        # 验证系统组分
        all_supported, rejected = self.validate_components(system)
        self.rejected_components.extend(rejected)
        
        # 检查是否为重复计算
        is_duplicate, previous_path, system_hash = self.check_duplicate_calculation(system)
        if is_duplicate:
            # 将重复计算信息添加到拒绝组分中
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'system',
                'component_name': system.name,
                'reason': f"重复计算，已有结果在: {previous_path}"
            })
            self.logger.warning(f"系统 {system.name} 已经计算过，跳过计算")
            return False
        
        # 准备目标目录
        destination_path = self.prepare_destination(system)
        
        # 切换到目标目录
        original_dir = os.getcwd()
        os.chdir(destination_path)
        
        try:
            # 准备溶剂分子
            self.prepare_solvents(system, destination_path)
            
            # 运行Packmol
            if not self.run_packmol(system, destination_path):
                return False
            
            # 准备Moltemplate文件
            if not self.prepare_moltemplate_files(system, destination_path):
                return False
            
            # 准备LAMMPS文件
            if not self.prepare_lammps_files(system, destination_path):
                return False
            
            # 提交作业
            if not self.submit_job(system, destination_path, submit_job):
                return False
            
            # 保存计算历史
            self.save_calculation_history(system_hash, destination_path)
            
            self.logger.info(f"Successfully completed simulation setup for {system.name}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error during simulation setup for {system.name}: {str(e)}")
            self.rejected_components.append({
                'system_name': system.name,
                'component_type': 'process',
                'component_name': 'simulation_setup',
                'reason': f"模拟设置过程中出错: {str(e)}"
            })
            return False
        finally:
            # 恢复原始目录
            os.chdir(original_dir)
            # 保存拒绝组分列表
            self.save_rejected_components() 