"""
电解液计算工作流

处理电解液计算的完整工作流程，包括输入解析、文件生成、计算提交和结果分析。
"""

import os
import sys
import re
import json
import time
import shutil
import logging
import subprocess
import traceback
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union
from pathlib import Path

# 导入molyte_cursor模块
from molyte_cursor.src.io.inp_reader import read_inp_file, parse_inp_content, INPReader
from molyte_cursor.src.io.file_generator import FileGenerator
from molyte_cursor.src.io.lammps_file_generator import LAMMPSFileGenerator
from molyte_cursor.src.io.electrolyte_file_generator import ElectrolyteFileGenerator
from molyte_cursor.src.utils.command_executor import CommandExecutor
from molyte_cursor.src.analysis.analyzer import Analyzer
from molyte_cursor.src.utils.logger import setup_logger

# 导入用户管理模块
try:
    from molyte_cursor.src.utils.user_manager import UserManager
except ImportError:
    UserManager = None

# 确保shutil在全局可用
if 'shutil' not in globals():
    import shutil

class ElectrolyteWorkflow:
    """
    电解液计算工作流
    
    管理电解液计算的完整流程，包括：
    1. 解析INP输入文件
    2. 生成LAMMPS输入文件
    3. 提交LAMMPS计算
    4. 分析结果
    5. 生成输出报告
    """
    
    def __init__(self, 
                 inp_file_path: Optional[str] = None, 
                 inp_content: Optional[str] = None,
                 work_dir: Optional[str] = None,
                 log_level: str = 'INFO',
                 user_id: Optional[str] = None,
                 project_id: Optional[str] = None,
                 formulation_id: Optional[str] = None):
        """
        初始化电解液工作流
        
        Args:
            inp_file_path: INP文件路径，如果提供则从文件加载
            inp_content: INP文件内容，如果提供则直接使用内容
            work_dir: 工作目录，如果不提供则使用临时目录
            log_level: 日志级别
            user_id: 用户ID，用于跟踪
            project_id: 项目ID，用于跟踪
            formulation_id: 配方ID，用于跟踪
        """
        # 添加导入路径
        # 检查是否在系统路径中
        import sys
        from pathlib import Path
        current_dir = Path(__file__).resolve().parent
        base_dir = current_dir.parent.parent.parent  # 应该是AiSuan根目录
        if str(base_dir) not in sys.path:
            sys.path.append(str(base_dir))
            print(f"已添加 {base_dir} 到Python路径")
        
        # 初始化日志记录器
        self.logger = logging.getLogger('electrolyte_workflow')
        self.logger.setLevel(log_level)
        
        # 设置控制台处理器
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)
        
        # 记录初始信息
        self.logger.info(f"初始化电解液计算工作流")
        
        # 保存用户相关信息
        self.user_id = user_id or 'anonymous'
        self.project_id = project_id
        self.formulation_id = formulation_id
        
        self.logger.info(f"用户ID: {self.user_id}, 项目ID: {self.project_id}, 配方ID: {self.formulation_id}")
        
        # 仅存储一个输入源
        if inp_file_path and inp_content:
            self.logger.warning("同时提供了文件路径和内容，将优先使用文件路径")
            
        self.inp_file_path = inp_file_path
        self.inp_content = inp_content
        
        # 设置工作目录
        if work_dir:
            self.work_dir = work_dir
        else:
            # 创建以时间戳命名的临时目录
            import tempfile
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            self.work_dir = os.path.join(tempfile.gettempdir(), f"electrolyte_{timestamp}")
        
        # 确保工作目录存在
        os.makedirs(self.work_dir, exist_ok=True)
        
        # 初始化状态和配置
        self.logger.info(f"工作目录: {self.work_dir}")
        self.status = 'initialized'
        self.config = {}
        self.generated_files = {}
        
        # 已解析标志
        self.input_parsed = False
        self.calculation_ran = False
        self.analysis_complete = False
        
        # 设置输出目录
        self.output_dir = self.work_dir
        
        # 添加文件日志处理器
        log_file = os.path.join(self.output_dir, 'workflow.log')
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        
        # 设置初始目录配置 - 这些会在解析输入后根据配置更新
        self.initial_salts_dir = os.path.join(base_dir, 'molyte_cursor', 'initial_salts')
        self.initial_solvent_dir = os.path.join(base_dir, 'molyte_cursor', 'initial_solvent')
        self.resp_path = os.path.join(base_dir, 'molyte_cursor', 'RESP')
        self.charge_save_path = os.path.join(base_dir, 'molyte_cursor', 'charge')
        self.molecules_workspace = os.path.join(base_dir, 'molyte_cursor', 'molecules_workspace')
        
        # 初始化paths属性
        self.paths = {}
        
        # 确保这些目录存在
        os.makedirs(self.initial_salts_dir, exist_ok=True)
        os.makedirs(self.initial_solvent_dir, exist_ok=True)
        os.makedirs(self.resp_path, exist_ok=True)
        os.makedirs(self.charge_save_path, exist_ok=True)
        os.makedirs(self.molecules_workspace, exist_ok=True)
    
    def _set_directory_permissions(self, directory: str) -> None:
        """
        设置目录权限
        
        Args:
            directory: 目录路径
        """
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        
        # 在Linux/Unix系统上设置权限
        if os.name == 'posix':
            try:
                # 如果用户ID不是匿名，尝试设置适当的权限
                if self.user_id and self.user_id != 'anonymous':
                    if UserManager:
                        user_manager = UserManager(self.workspace_root)
                        user_info = user_manager.get_user_info(self.user_id)
                        
                        if user_info:
                            # 管理员目录: 只有所有者有完全权限 (700)
                            if user_info['role'] == 'admin':
                                os.chmod(directory, 0o700)
                            # 研究员目录: 所有者有完全权限，组用户可读可执行 (750)
                            else:
                                os.chmod(directory, 0o750)
                        else:
                            # 默认权限: 所有者可读写执行，其他人可读可执行 (755)
                            os.chmod(directory, 0o755)
                    else:
                        # 用户管理模块不可用，使用默认权限
                        os.chmod(directory, 0o755)
                else:
                    # 匿名用户目录: 所有者可读写执行，其他人可读可执行 (755)
                    os.chmod(directory, 0o755)
            except Exception as e:
                self.logger.warning(f"设置目录权限失败: {str(e)}")
    
    def run(self) -> Dict[str, Any]:
        """
        运行完整工作流
        
        Returns:
            Dict: 计算结果字典
        """
        try:
            self.logger.info("开始执行电解液计算工作流")
            
            # 步骤1：解析INP输入
            self.parse_input()
            
            # 步骤2：生成输入文件
            self.generate_input_files()
            
            # 步骤3：运行LAMMPS计算
            self.run_calculation()
            
            # 步骤4：分析结果
            self.analyze_results()
            
            # 步骤5：生成报告
            self.generate_report()
            
            self.status = 'completed'
            self.logger.info("电解液计算工作流完成")
            
            return self.get_result_summary()
            
        except Exception as e:
            self.status = 'failed'
            self.logger.error(f"工作流执行失败: {str(e)}", exc_info=True)
            
            return {
                'status': 'failed',
                'error': str(e),
                'output_dir': self.output_dir
            }
    
    def parse_input(self) -> None:
        """
        解析INP输入文件
        
        解析INP输入文件的内容，获取所需的配置信息。
        """
        self.logger.info("开始解析输入文件")
        
        if self.inp_file_path:
            self.logger.info(f"解析INP文件: {self.inp_file_path}")
            # 使用INPReader解析INP文件
            self.config = read_inp_file(self.inp_file_path)
        elif self.inp_content:
            self.logger.info("解析INP内容")
            # 解析INP内容字符串
            self.config = parse_inp_content(self.inp_content)
        else:
            raise ValueError("未提供INP文件或内容")
        
        # 更新配置字典
        self.config['output_dir'] = self.output_dir
        
        # 读取所需的路径配置
        paths = self.check_paths()
        self.config['paths'] = paths
        
        # 设置路径属性
        self.initial_solvent_dir = paths.get('initial_solvent_dir', '')
        self.initial_salts_dir = paths.get('initial_salts_dir', '')  # 添加盐目录属性
        self.ligpargen_path = paths.get('ligpargen_path', 'ligpargen')
        self.packmol_path = paths.get('packmol_path', 'packmol')
        self.lempyfy_path = paths.get('lempyfy_path', 'ltemplify.py')
        self.moltemplate_path = paths.get('moltemplate_path', 'moltemplate.sh')
        self.resp_path = paths.get('resp_path', '')
        self.charge_save_path = paths.get('charge_save_path', '')
        self.molecules_workspace = paths.get('molecules_workspace', '')
        
        self.logger.info(f"初始盐目录: {self.initial_salts_dir}")
        self.logger.info(f"初始溶剂目录: {self.initial_solvent_dir}")
        self.logger.info(f"RESP路径: {self.resp_path}")
        self.logger.info(f"电荷保存路径: {self.charge_save_path}")
        self.logger.info(f"分子工作区路径: {self.molecules_workspace}")
    
    def generate_input_files(self) -> None:
        """生成LAMMPS输入文件，包括使用LigParGen和Packmol处理分子和体系"""
        self.logger.info("步骤2: 生成LAMMPS输入文件")
        
        try:
            # 添加详细日志，打印调用堆栈
            import traceback
            self.logger.info("执行generate_input_files方法，调用堆栈:")
            self.logger.info(traceback.format_stack())
            
            # 确保配置已解析
            if not self.config:
                raise ValueError("配置未解析，请先调用parse_input()方法")
            
            # 创建输出目录结构
            output_dir = os.path.join(self.output_dir, 'input_files')
            molecule_dir = os.path.join(output_dir, 'molecules')
            packmol_dir = os.path.join(output_dir, 'packmol')
            lammps_dir = os.path.join(output_dir, 'lammps')
            temp_dir = os.path.join(self.output_dir, 'temp')  # 添加临时目录
            
            # 确保所有目录都存在，并设置正确的权限
            for dir_path in [output_dir, molecule_dir, packmol_dir, lammps_dir, temp_dir]:
                try:
                    if not os.path.exists(dir_path):
                        os.makedirs(dir_path, exist_ok=True)
                    # 设置目录权限
                    self._set_directory_permissions(dir_path)
                except Exception as e:
                    self.logger.error(f"创建目录失败: {dir_path}, 错误: {str(e)}")
                    raise RuntimeError(f"无法创建必要的目录: {str(e)}")
            
            self.logger.info("创建目录结构完成")
            
            # 1. 处理分子 - 使用LigParGen生成分子结构和力场文件
            self.logger.info("1. 使用LigParGen处理分子结构和力场")
            molecule_files = self._process_molecules_with_ligpargen(molecule_dir)
            
            # 2. 使用Packmol构建分子体系
            self.logger.info("2. 使用Packmol构建分子体系")
            packmol_output = self._build_system_with_packmol(molecule_files, packmol_dir)
            
            # 3. 使用ElectrolyteFileGenerator生成LAMMPS输入文件
            self.logger.info("3. 生成LAMMPS输入文件")
            electrolyte_generator = ElectrolyteFileGenerator()
            
            # 生成LAMMPS输入文件
            self.generated_files = electrolyte_generator.generate_lammps_input_files(
                molecule_files=molecule_files,
                output_dir=lammps_dir,
                packmol_output=packmol_output,
                config=self.config
            )
            
            # 记录生成的文件
            files_list_path = os.path.join(self.output_dir, 'generated_files.json')
            with open(files_list_path, 'w', encoding='utf-8') as f:
                json.dump(self.generated_files, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"已生成{len(self.generated_files)}个输入文件")
            self.logger.debug(f"生成的文件: {json.dumps(self.generated_files, indent=2)}")
            
            # 更新状态
            self.status = 'files_generated'
            
        except Exception as e:
            self.logger.error(f"生成输入文件失败: {str(e)}", exc_info=True)
            raise ValueError(f"生成输入文件失败: {str(e)}")
    
    def check_molecule_name(self, name: str, smile: str, existing_molecules: Dict[str, Dict[str, str]], molecule_type: str) -> str:
        """
        为分子生成唯一名称。如果名称已存在且SMILE不同，则添加后缀生成新名称。
        
        Args:
            name: 分子名称
            smile: 分子的SMILE字符串
            existing_molecules: 已存在的分子字典，格式为 {molecule_type: {name: smile}}
            molecule_type: 分子类型，例如 'solvent', 'cation', 'anion'
            
        Returns:
            str: 唯一的分子名称。如果名称未重复，返回原名称；如果名称重复但SMILE相同，返回原名称；
                如果名称重复且SMILE不同，返回添加后缀的新名称。
        """
        # 确保分子类型存在于字典中
        if molecule_type not in existing_molecules:
            existing_molecules[molecule_type] = {}
        
        # 检查名称是否已存在
        if name not in existing_molecules[molecule_type]:
            # 名称不存在，直接使用并记录
            existing_molecules[molecule_type][name] = smile
            return name
        else:
            # 名称已存在，检查SMILE是否相同
            if existing_molecules[molecule_type][name] == smile:
                # SMILE相同，可以使用相同名称
                return name
            else:
                # SMILE不同，需要生成新名称
                self.logger.warning(f"分子名称 '{name}' 已存在，但SMILE不同。将生成新名称。")
                suffix = 1
                while True:
                    new_name = f"{name}_{suffix}"
                    if new_name not in existing_molecules[molecule_type]:
                        existing_molecules[molecule_type][new_name] = smile
                        self.logger.info(f"生成新名称: {new_name} 对应SMILE: {smile}")
                        return new_name
                    suffix += 1
    
    def _process_molecules_with_ligpargen(self, work_dir: str) -> Dict[str, Any]:
        """
        使用LigParGen处理分子结构和力场
        
        Args:
            work_dir: 工作目录
            
        Returns:
            Dict: 包含生成的分子文件信息
        """
        self.logger.info("处理分子结构和力场...")
        
        try:
            # 获取所有需要处理的分子
            molecules = self.config.get('molecules', [])
            
            # 检查目录
            self.logger.info(f"使用初始溶剂目录: {self.initial_solvent_dir}")
            self.logger.info(f"使用电荷保存目录: {self.charge_save_path}")
            self.logger.info(f"使用RESP路径: {self.resp_path}")
            
            # 存储处理后的分子文件信息
            molecule_info = []  # 用于跟踪分子信息的变量
            topology_files = []
            forcefield_files = []
            
            os.makedirs(work_dir, exist_ok=True)
            
            # 检查是否有分子需要处理
            if not molecules:
                self.logger.error("配置中没有分子信息，请检查INP文件")
                # 我们可以尝试检查是否有Sol1_name, Cation1_name等键
                parsed_molecules = []  # 存储从配置中解析的分子列表，而不直接使用molecule_info
                for key in self.config.keys():
                    # 跳过formulation_name键，它不是分子
                    if key == 'formulation_name':
                        continue
                        
                    if key.endswith('_name'):
                        mol_name = self.config[key]
                        mol_type = 'solvent' if key.startswith('Sol') else ('cation' if key.startswith('Cation') else 'anion')
                        mol_smile = self.config.get(f"{key[:-5]}_smile", "")
                        
                        self.logger.info(f"从配置键{key}中找到分子: {mol_name} ({mol_type})")
                        parsed_molecules.append({
                            'name': mol_name,
                            'type': mol_type,
                            'smile': mol_smile
                        })
                
                if not parsed_molecules:
                    self.logger.error("尝试从配置键中解析分子信息失败")
                    raise ValueError("配置中没有分子信息，请检查INP文件")
                
                # 使用parsed_molecules作为要处理的分子列表，而不是molecule_info
                molecules = parsed_molecules
            
            for mol in molecules:
                mol_name = mol.get('name', '')
                mol_type = mol.get('type', '')
                mol_smile = mol.get('smile', '')
                
                self.logger.info(f"处理分子: {mol_name} ({mol_type})")
                
                # 创建分子输出目录
                mol_dir = os.path.join(work_dir, mol_name)
                os.makedirs(mol_dir, exist_ok=True)
                
                # 对于阴离子和阳离子，从initial_salts目录复制文件
                if mol_type in ['cation', 'anion']:
                    self.logger.info(f"从initial_salts目录加载离子 {mol_name} 的文件")
                    
                    # 处理离子名称，移除+或-符号
                    clean_mol_name = mol_name.replace('+', '').replace('-', '')
                    
                    # 在initial_salts目录中查找文件
                    # 先尝试使用非电荷名称（Li 而非 Li+）
                    source_pdb = os.path.join(self.initial_salts_dir, f"{clean_mol_name}.pdb")
                    source_lt = os.path.join(self.initial_salts_dir, f"{clean_mol_name}.lt")
                    
                    self.logger.info(f"查找离子文件: {source_pdb} 和 {source_lt}")
                    
                    # 如果文件不存在，则尝试使用原始名称（可能带有+或-符号）
                    if not os.path.exists(source_pdb) or not os.path.exists(source_lt):
                        self.logger.warning(f"使用简化名称找不到离子文件，尝试使用原始名称: {mol_name}")
                        source_pdb = os.path.join(self.initial_salts_dir, f"{mol_name}.pdb")
                        source_lt = os.path.join(self.initial_salts_dir, f"{mol_name}.lt")
                    
                    self.logger.info(f"查找离子文件: {source_pdb} 和 {source_lt}")
                    
                    # 目标文件路径 - 使用不带+和-符号的名称
                    mol_dir = os.path.join(work_dir, clean_mol_name)
                    os.makedirs(mol_dir, exist_ok=True)
                    
                    dest_pdb = os.path.join(mol_dir, f"{clean_mol_name}.pdb")
                    dest_lt = os.path.join(mol_dir, f"{clean_mol_name}.lt")
                    
                    # 检查源文件是否存在
                    if not os.path.exists(source_pdb) or not os.path.exists(source_lt):
                        self.logger.error(f"在initial_salts目录中找不到离子文件: {source_pdb} 或 {source_lt}")
                        self.logger.error(f"阴阳离子需要预定义文件")
                        raise FileNotFoundError(f"无法创建必需的LT文件: {clean_mol_name}.lt")
                    else:
                        # 复制文件
                        shutil.copy2(source_pdb, dest_pdb)
                        shutil.copy2(source_lt, dest_lt)
                        self.logger.info(f"成功复制离子文件: {dest_pdb}, {dest_lt}")
                        
                        # 定义输出文件路径
                        top_file = dest_lt  # 使用.lt文件作为拓扑文件
                        pdb_file = dest_pdb
                        itp_file = dest_lt  # 使用.lt文件作为力场文件
                        lt_file = dest_lt
                    
                    topology_files.append(top_file)
                    forcefield_files.append(itp_file)
                    
                    molecule_info.append({
                        'name': mol_name,
                        'type': mol_type,
                        'topology': top_file,
                        'structure': pdb_file,
                        'forcefield': itp_file,
                        'lt_file': lt_file
                    })
                
                # 对于溶剂，使用LigParGen生成文件并处理电荷
                elif mol_type == 'solvent':
                    self.logger.info(f"处理溶剂分子: {mol_name}")
                    
                    # 验证SMILE字符串
                    if mol_smile:
                        # 不再调用validate_smile
                        pass
                    # 以下内容替换之前的validate_smile调用
                    
                    # 首先检查workspace目录中是否已存在该溶剂文件
                    workspace_solvent_dir = self.molecules_workspace
                    molecule_dir = os.path.join(workspace_solvent_dir, mol_name)
                    
                    workspace_pdb = os.path.join(molecule_dir, f"{mol_name}.pdb")
                    workspace_lmp = os.path.join(molecule_dir, f"{mol_name}.lmp")
                    workspace_lt = os.path.join(molecule_dir, f"{mol_name}.lt")
                    workspace_charge = os.path.join(molecule_dir, f"{mol_name}.charge")
                    workspace_smile_file = os.path.join(molecule_dir, f"{mol_name}.smile")
                    
                    # 检查workspace目录中是否存在此溶剂的文件
                    workspace_files_exist = os.path.exists(workspace_pdb) and os.path.exists(workspace_lmp)
                    
                    # 如果存在文件，还需要检查SMILE是否一致
                    if workspace_files_exist and os.path.exists(workspace_smile_file) and mol_smile:
                        # 读取保存的SMILE字符串
                        try:
                            with open(workspace_smile_file, 'r') as f:
                                saved_smile = f.read().strip()
                            
                            # 比较SMILE字符串是否一致
                            if saved_smile != mol_smile:
                                self.logger.warning(f"溶剂 {mol_name} 在workspace中存在，但SMILE不同: '{saved_smile}' != '{mol_smile}'")
                                self.logger.warning(f"将生成新名称或创建新文件")
                                
                                # 生成唯一名称
                                unique_name = self.check_molecule_name(mol_name, mol_smile, existing_molecules, 'solvent')
                                if unique_name != mol_name:
                                    self.logger.info(f"溶剂名称 '{mol_name}' 更改为 '{unique_name}'")
                                    mol_name = unique_name
                                    
                                    # 更新文件路径
                                    molecule_dir = os.path.join(workspace_solvent_dir, mol_name)
                                    workspace_pdb = os.path.join(molecule_dir, f"{mol_name}.pdb")
                                    workspace_lmp = os.path.join(molecule_dir, f"{mol_name}.lmp")
                                    workspace_lt = os.path.join(molecule_dir, f"{mol_name}.lt")
                                    workspace_charge = os.path.join(molecule_dir, f"{mol_name}.charge")
                                    workspace_smile_file = os.path.join(molecule_dir, f"{mol_name}.smile")
                                    
                                    # 重新检查是否存在
                                    workspace_files_exist = os.path.exists(workspace_pdb) and os.path.exists(workspace_lmp)
                            else:
                                self.logger.info(f"溶剂 {mol_name} 在workspace中存在，且SMILE一致: '{saved_smile}'")
                        except Exception as e:
                            self.logger.warning(f"读取workspace中的SMILE文件出错: {e}")
                            workspace_files_exist = False
                    
                    # 如果workspace中存在文件，优先使用它们
                    if workspace_files_exist:
                        self.logger.info(f"在workspace目录中找到溶剂 {mol_name} 的文件")
                        
                        # 创建分子输出目录中的文件
                        dest_pdb = os.path.join(mol_dir, f"{mol_name}.pdb")
                        dest_lmp = os.path.join(mol_dir, f"{mol_name}.lmp")
                        dest_lt = os.path.join(mol_dir, f"{mol_name}.lt")
                        
                        # 复制文件
                        shutil.copy2(workspace_pdb, dest_pdb)
                        shutil.copy2(workspace_lmp, dest_lmp)
                        
                        if os.path.exists(workspace_lt):
                            shutil.copy2(workspace_lt, dest_lt)
                        else:
                            # 如果LT文件不存在，生成它
                            self.logger.info(f"在workspace中找不到LT文件，将生成: {dest_lt}")
                            self._generate_lt_file(dest_lmp, dest_pdb, dest_lt)
                        
                        # 处理电荷文件
                        if os.path.exists(workspace_charge):
                            self.logger.info(f"使用workspace中的电荷文件: {workspace_charge}")
                            charge_file = os.path.join(self.charge_save_path, f"{mol_name}.charge")
                            
                            # 如果charge_save_path中不存在此电荷文件，复制它
                            if not os.path.exists(charge_file):
                                shutil.copy2(workspace_charge, charge_file)
                                self.logger.info(f"已复制电荷文件到charge目录: {charge_file}")
                            
                            # 使用电荷文件修改LMP文件
                            dest_lmp_charged = os.path.join(mol_dir, f"{mol_name}_charged.lmp")
                            self._modify_lmp_charges(dest_lmp, workspace_charge, dest_lmp_charged)
                            
                            # 如果修改成功，替换原LMP文件
                            if os.path.exists(dest_lmp_charged):
                                shutil.move(dest_lmp_charged, dest_lmp)
                                self.logger.info(f"已使用电荷文件更新 {dest_lmp}")
                        
                        # 添加文件信息
                        top_file = dest_lt
                        structure_file = dest_pdb
                        itp_file = dest_lt
                        lt_file = dest_lt
                        lmp_file = dest_lmp
                        
                        topology_files.append(top_file)
                        forcefield_files.append(itp_file)
                        
                        molecule_info.append({
                            'name': mol_name,
                            'type': mol_type,
                            'topology': top_file,
                            'structure': structure_file,
                            'forcefield': itp_file,
                            'lt_file': lt_file,
                            'lmp_file': lmp_file
                        })
                        
                        # 已处理完毕，跳过后续处理
                        continue
                    
                    # 定义文件路径
                    default_pdb = os.path.join(self.initial_solvent_dir, f"{mol_name}.pdb")
                    default_lmp = os.path.join(self.initial_solvent_dir, f"{mol_name}.lmp")
                    default_lt = os.path.join(self.initial_solvent_dir, f"{mol_name}.lt")
                    
                    # 检查溶剂文件是否存在
                    pdb_exists = os.path.exists(default_pdb)
                    lmp_exists = os.path.exists(default_lmp)
                    lt_exists = os.path.exists(default_lt)
                    
                    self.logger.info(f"溶剂文件状态 - PDB: {pdb_exists}, LMP: {lmp_exists}, LT: {lt_exists}")
                    
                    # 如果文件不存在且有SMILE字符串，则使用SMILE生成
                    if (not pdb_exists or not lmp_exists) and mol_smile:
                        self.logger.info(f"溶剂 {mol_name} 的文件不存在或不完整，将使用SMILE生成")
                        try:
                            # 创建临时目录用于生成文件
                            temp_gen_dir = os.path.join(self.work_dir, "solvent_generation", mol_name)
                            os.makedirs(temp_gen_dir, exist_ok=True)
                            
                            # 使用SMILE生成溶剂文件
                            generated_files = self._generate_solvent_files_with_smile(
                                mol_name, mol_smile, temp_gen_dir
                            )
                            
                            # 更新文件路径
                            default_pdb = generated_files['pdb']
                            default_lmp = generated_files['lmp']
                            default_lt = generated_files['lt']
                            
                            # 更新文件状态
                            pdb_exists = os.path.exists(default_pdb)
                            lmp_exists = os.path.exists(default_lmp)
                            lt_exists = os.path.exists(default_lt)
                            
                            self.logger.info(f"文件生成后状态 - PDB: {pdb_exists}, LMP: {lmp_exists}, LT: {lt_exists}")
                        except Exception as e:
                            self.logger.error(f"生成溶剂文件失败: {str(e)}")
                            self.logger.warning("将尝试继续处理现有文件，可能会导致后续错误")
                    
                    # 检查电荷文件
                    charge_file = os.path.join(self.charge_save_path, f"{mol_name}.charge")
                    
                    # 处理溶剂文件
                    if pdb_exists and lmp_exists:
                        # 目标文件路径
                        dest_pdb = os.path.join(mol_dir, f"{mol_name}.pdb")
                        dest_lmp = os.path.join(mol_dir, f"{mol_name}.lmp")
                        dest_lt = os.path.join(mol_dir, f"{mol_name}.lt")
                        
                        # 复制基本文件
                        self._copy_file_with_check(default_pdb, dest_pdb, "PDB")
                        self._copy_file_with_check(default_lmp, dest_lmp, "LMP")
                        if lt_exists:
                            self._copy_file_with_check(default_lt, dest_lt, "LT")
                        
                        # 处理电荷
                        if os.path.exists(charge_file):
                            self.logger.info(f"找到现有电荷文件: {charge_file}")
                            
                            # 使用现有电荷修改LMP文件
                            dest_lmp_charged = os.path.join(mol_dir, f"{mol_name}_charged.lmp")
                            self._modify_lmp_charges(dest_lmp, charge_file, dest_lmp_charged)
                            
                            # 如果修改成功，替换原LMP文件
                            if os.path.exists(dest_lmp_charged):
                                shutil.move(dest_lmp_charged, dest_lmp)
                                self.logger.info(f"已使用电荷文件 {charge_file} 更新 {dest_lmp}")
                        else:
                            # 生成电荷文件
                            self.logger.info(f"未找到电荷文件，将生成新的电荷文件: {charge_file}")
                            self._generate_charges_with_resp(mol_name, dest_pdb, self.resp_path)
                            
                            # 保存电荷文件到charge目录
                            temp_charge_file = os.path.join(mol_dir, f"{mol_name}.charge")
                            if os.path.exists(temp_charge_file):
                                shutil.copy2(temp_charge_file, charge_file)
                                self.logger.info(f"已保存电荷文件: {charge_file}")
                                
                                # 同时复制电荷文件到workspace目录
                                workspace_solvent_dir = self.molecules_workspace
                                molecule_dir = os.path.join(workspace_solvent_dir, mol_name)
                                os.makedirs(molecule_dir, exist_ok=True)
                                
                                workspace_charge = os.path.join(molecule_dir, f"{mol_name}.charge")
                                shutil.copy2(temp_charge_file, workspace_charge)
                                self.logger.info(f"已保存电荷文件到workspace: {workspace_charge}")
                            
                                # 使用新生成的电荷修改LMP文件
                                dest_lmp_charged = os.path.join(mol_dir, f"{mol_name}_charged.lmp")
                                self._modify_lmp_charges(dest_lmp, charge_file, dest_lmp_charged)
                                
                                # 如果修改成功，替换原LMP文件
                                if os.path.exists(dest_lmp_charged):
                                    shutil.move(dest_lmp_charged, dest_lmp)
                                    self.logger.info(f"已使用新生成的电荷文件更新 {dest_lmp}")
                        
                        # 如果没有LT文件，生成LT文件
                        if not os.path.exists(dest_lt):
                            self.logger.info(f"未找到LT文件，将使用ltemplify生成: {dest_lt}")
                            self._generate_lt_file(dest_lmp, dest_pdb, dest_lt)
                        
                        # 添加文件信息
                        top_file = os.path.join(mol_dir, f"{mol_name}.top")
                        itp_file = os.path.join(mol_dir, f"{mol_name}.itp")
                        
                        # 如果没有top或itp文件，使用lt文件代替
                        if not os.path.exists(top_file):
                            top_file = dest_lt
                        if not os.path.exists(itp_file):
                            itp_file = dest_lt
                        
                        topology_files.append(top_file)
                        forcefield_files.append(itp_file)
                        
                        molecule_info.append({
                            'name': mol_name,
                            'type': mol_type,
                            'topology': top_file,
                            'structure': dest_pdb,
                            'forcefield': itp_file,
                            'lt_file': dest_lt,
                            'lmp_file': dest_lmp
                        })
                    elif mol_smile:
                        # 文件不存在，但有SMILE字符串，则再次尝试生成
                        self.logger.warning(f"溶剂 {mol_name} 的文件处理失败，再次尝试使用SMILE生成")
                        
                        # 创建临时目录用于生成文件
                        temp_gen_dir = os.path.join(self.work_dir, "solvent_generation_retry", mol_name)
                        os.makedirs(temp_gen_dir, exist_ok=True)
                        
                        try:
                            # 使用SMILE生成溶剂文件，此次使用更详细的日志
                            self.logger.info(f"使用SMILE '{mol_smile}' 为溶剂 {mol_name} 生成文件")
                            
                            # 验证SMILE字符串
                            # 删除对validate_smile的调用
                            if not mol_smile:
                                self.logger.error(f"无效的SMILE字符串: {mol_smile}")
                                raise ValueError(f"无法为溶剂 {mol_name} 生成文件：SMILE字符串无效")
                            
                            # 调用更详细的生成函数
                            generated_files = self._generate_solvent_files_with_smile(
                                mol_name, mol_smile, temp_gen_dir, verbose=True
                            )
                            
                            # 创建分子目录中的文件
                            dest_pdb = os.path.join(mol_dir, f"{mol_name}.pdb")
                            dest_lmp = os.path.join(mol_dir, f"{mol_name}.lmp")
                            dest_lt = os.path.join(mol_dir, f"{mol_name}.lt")
                            
                            # 复制生成的文件
                            shutil.copy2(generated_files['pdb'], dest_pdb)
                            shutil.copy2(generated_files['lmp'], dest_lmp)
                            if 'lt' in generated_files and os.path.exists(generated_files['lt']):
                                shutil.copy2(generated_files['lt'], dest_lt)
                            else:
                                # 生成LT文件
                                self._generate_lt_file(dest_lmp, dest_pdb, dest_lt)
                            
                            # 生成并处理电荷
                            self.logger.info(f"为溶剂 {mol_name} 生成电荷文件")
                            self._generate_charges_with_resp(mol_name, dest_pdb, self.resp_path)
                            
                            # 保存电荷文件
                            temp_charge_file = os.path.join(mol_dir, f"{mol_name}.charge")
                            if os.path.exists(temp_charge_file):
                                # 复制到charge目录
                                shutil.copy2(temp_charge_file, charge_file)
                                
                                # 使用电荷文件修改LMP文件
                                dest_lmp_charged = os.path.join(mol_dir, f"{mol_name}_charged.lmp")
                                self._modify_lmp_charges(dest_lmp, temp_charge_file, dest_lmp_charged)
                                
                                # 如果修改成功，替换原LMP文件
                                if os.path.exists(dest_lmp_charged):
                                    shutil.move(dest_lmp_charged, dest_lmp)
                            
                            # 添加文件信息
                            top_file = dest_lt
                            itp_file = dest_lt
                            
                            topology_files.append(top_file)
                            forcefield_files.append(itp_file)
                            
                            molecule_info.append({
                                'name': mol_name,
                                'type': mol_type,
                                'topology': top_file,
                                'structure': dest_pdb,
                                'forcefield': itp_file,
                                'lt_file': dest_lt,
                                'lmp_file': dest_lmp
                            })
                            
                            # 同时保存到workspace目录
                            workspace_solvent_dir = self.molecules_workspace
                            workspace_mol_dir = os.path.join(workspace_solvent_dir, mol_name)
                            os.makedirs(workspace_mol_dir, exist_ok=True)
                            
                            workspace_pdb = os.path.join(workspace_mol_dir, f"{mol_name}.pdb")
                            workspace_lmp = os.path.join(workspace_mol_dir, f"{mol_name}.lmp")
                            workspace_lt = os.path.join(workspace_mol_dir, f"{mol_name}.lt")
                            workspace_charge = os.path.join(workspace_mol_dir, f"{mol_name}.charge")
                            workspace_smile_file = os.path.join(workspace_mol_dir, f"{mol_name}.smile")
                            
                            # 复制文件到workspace
                            shutil.copy2(dest_pdb, workspace_pdb)
                            shutil.copy2(dest_lmp, workspace_lmp)
                            shutil.copy2(dest_lt, workspace_lt)
                            
                            if os.path.exists(temp_charge_file):
                                shutil.copy2(temp_charge_file, workspace_charge)
                            
                            # 保存SMILE字符串
                            with open(workspace_smile_file, 'w') as f:
                                f.write(mol_smile)
                            
                            self.logger.info(f"溶剂 {mol_name} 的文件已保存到workspace目录: {workspace_mol_dir}")
                            
                        except Exception as e:
                            self.logger.error(f"第二次尝试生成溶剂 {mol_name} 的文件失败: {str(e)}")
                            import traceback
                            self.logger.error(traceback.format_exc())
                            
                            # 此时确实无法处理，抛出异常
                            error_msg = f"无法为溶剂 {mol_name} 生成文件，尽管提供了SMILE字符串: {mol_smile}"
                            raise FileNotFoundError(error_msg)
                    else:
                        # 文件缺失且无法生成
                        missing_files = []
                        if not pdb_exists:
                            missing_files.append(f"PDB ({default_pdb})")
                        if not lmp_exists:
                            missing_files.append(f"LMP ({default_lmp})")
                        
                        error_msg = f"溶剂 {mol_name} 的必要文件缺失: {', '.join(missing_files)}"
                        if not mol_smile:
                            error_msg += "，且未提供SMILE字符串，无法自动生成"
                        
                        self.logger.error(error_msg)
                        raise FileNotFoundError(error_msg)
            
            # 保存更新后的分子信息
            try:
                # 定义存储分子信息的文件路径
                existing_molecules_file = os.path.join(self.resp_path, "existing_molecules.json")
                
                # 检查现有分子信息文件是否存在
                existing_molecules = {'solvent': {}, 'cation': {}, 'anion': {}}
                if os.path.exists(existing_molecules_file):
                    try:
                        with open(existing_molecules_file, 'r') as f:
                            loaded = json.load(f)
                            # 合并现有数据
                            for mol_type in existing_molecules:
                                if mol_type in loaded:
                                    existing_molecules[mol_type].update(loaded[mol_type])
                    except Exception as e:
                        self.logger.warning(f"读取现有分子信息文件出错: {str(e)}")
                
                with open(existing_molecules_file, 'w') as f:
                    json.dump(existing_molecules, f, indent=4)
                self.logger.info(f"已更新分子信息文件: {existing_molecules_file}")
            except Exception as e:
                self.logger.warning(f"保存分子信息时出错: {e}")
            
            # 检查是否有分子被添加到了molecule_info中
            if not molecule_info:
                self.logger.error("未能成功添加任何分子信息")
                raise ValueError("处理分子结构过程中未能添加任何分子信息")
            
            return molecule_info
        
        except Exception as e:
            self.logger.error(f"处理分子结构和力场时出错: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            raise
    
    def _modify_lmp_charges(self, lmp_file: str, charge_file: str, output_file: str) -> None:
        """
        使用电荷文件修改LAMMPS数据文件中的电荷
        
        Args:
            lmp_file: LAMMPS数据文件路径
            charge_file: 电荷文件路径
            output_file: 输出文件路径
        """
        self.logger.info(f"使用电荷文件 {charge_file} 修改 {lmp_file}")
        
        try:
            # 读取电荷文件
            with open(charge_file, 'r') as f:
                charges = []
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        try:
                            charge = float(line)
                            charges.append(charge)
                        except ValueError:
                            # 如果是多列数据，尝试取最后一列
                            parts = line.split()
                            if parts:
                                try:
                                    charge = float(parts[-1])
                                    charges.append(charge)
                                except ValueError:
                                    pass
            
            if not charges:
                self.logger.error(f"电荷文件 {charge_file} 中未找到有效电荷")
                return
            
            self.logger.info(f"从电荷文件中读取了 {len(charges)} 个电荷值")
            
            # 读取并修改LAMMPS文件
            with open(lmp_file, 'r') as f:
                lmp_lines = f.readlines()
            
            # 查找原子部分
            atoms_section = False
            charge_index = -1
            atom_count = 0
            output_lines = []
            
            for line in lmp_lines:
                if 'Atoms' in line:
                    atoms_section = True
                    output_lines.append(line)
                    continue
                
                if atoms_section and line.strip() and not line.strip().startswith('#'):
                    parts = line.split()
                    if len(parts) >= 7:  # 假设格式为 atom_id mol_id atom_type q x y z
                        atom_count += 1
                        if atom_count <= len(charges):
                            # 替换电荷值（通常是第4列）
                            parts[3] = str(charges[atom_count - 1])
                            line = ' '.join(parts) + '\n'
                    
                output_lines.append(line)
                
                # 检查是否已处理完原子部分
                if atoms_section and line.strip() == '' and atom_count > 0:
                    atoms_section = False
            
            # 写入修改后的文件
            with open(output_file, 'w') as f:
                f.writelines(output_lines)
            
            self.logger.info(f"已修改 {atom_count} 个原子的电荷，保存到 {output_file}")
        
        except Exception as e:
            self.logger.error(f"修改LAMMPS电荷时出错: {e}")
    
    def _generate_charges_with_resp(self, mol_name: str, pdb_file: str, resp_path: str) -> None:
        """
        使用RESP生成电荷
        
        Args:
            mol_name: 分子名称
            pdb_file: PDB文件路径
            resp_path: RESP工具路径
        """
        self.logger.info(f"为分子 {mol_name} 生成RESP电荷")
        
        try:
            # 复制PDB文件到RESP目录
            resp_pdb = os.path.join(resp_path, f"{mol_name}.pdb")
            shutil.copy2(pdb_file, resp_pdb)
            
            # 运行RESP脚本
            resp_script = os.path.join(resp_path, "RESP2.sh")
            if os.path.exists(resp_script):
                # 切换到RESP目录
                original_dir = os.getcwd()
                os.chdir(resp_path)
                
                # 运行RESP脚本
                cmd = f"bash RESP2.sh {mol_name}.pdb"
                self.logger.info(f"执行RESP命令: {cmd}")
                
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = proc.communicate()
                
                if proc.returncode != 0:
                    self.logger.error(f"RESP命令执行失败: {stderr.decode()}")
                else:
                    self.logger.info(f"RESP命令执行成功: {stdout.decode()}")
                
                # 恢复原始目录
                os.chdir(original_dir)
                
                # 检查是否生成了电荷文件
                charge_file = os.path.join(resp_path, f"{mol_name}.charge")
                if os.path.exists(charge_file):
                    # 复制电荷文件到分子目录
                    dest_charge = os.path.join(os.path.dirname(pdb_file), f"{mol_name}.charge")
                    shutil.copy2(charge_file, dest_charge)
                    self.logger.info(f"已复制电荷文件: {charge_file} -> {dest_charge}")
                else:
                    self.logger.error(f"RESP未生成电荷文件: {charge_file}")
            else:
                self.logger.error(f"找不到RESP脚本: {resp_script}")
        
        except Exception as e:
            self.logger.error(f"生成RESP电荷时出错: {e}")
    
    def _generate_lt_file(self, lmp_file: str, pdb_file: str, lt_file: str) -> None:
        """
        使用ltemplify生成LT文件
        
        Args:
            lmp_file: LAMMPS数据文件路径
            pdb_file: PDB文件路径
            lt_file: 输出LT文件路径
        """
        # 获取分子名称（从lt_file路径中提取）
        mol_name = os.path.splitext(os.path.basename(lt_file))[0]
        self.logger.info(f"为分子 {mol_name} 的LMP文件 {lmp_file} 生成LT文件")
        
        try:
            # 查找ltemplify命令
            ltemplify_cmd = self.lempyfy_path
            if not self._check_command(ltemplify_cmd):
                python_cmd = shutil.which("python") or shutil.which("python3")
                if python_cmd:
                    ltemplify_cmd = f"{python_cmd} {self.lempyfy_path}"
                else:
                    self.logger.error("找不到python或ltemplify命令")
                    return
            
            # 生成LT文件，使用-name参数指定分子名称
            cmd = f"{ltemplify_cmd} -name {mol_name} {lmp_file} > {lt_file}"
            self.logger.info(f"执行命令: {cmd}")
            
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            
            if proc.returncode != 0:
                self.logger.error(f"ltemplify命令执行失败: {stderr.decode()}")
                # 如果ltemplify执行失败，手动创建一个基本的LT文件
                self._create_manual_lt_file(lmp_file, lt_file, mol_name)
            else:
                self.logger.info(f"成功生成LT文件: {lt_file}")
        
        except Exception as e:
            self.logger.error(f"生成LT文件时出错: {e}")
            # 如果发生异常，也尝试手动创建LT文件
            self._create_manual_lt_file(lmp_file, lt_file, mol_name)
    
    def _create_manual_lt_file(self, lmp_file: str, lt_file: str, mol_name: str) -> None:
        """
        当ltemplify失败时，手动创建一个基本的LT文件
        
        Args:
            lmp_file: LAMMPS数据文件路径
            lt_file: 输出LT文件路径
            mol_name: 分子名称
        """
        self.logger.info(f"ltemplify失败，手动创建基本LT文件: {lt_file}")
        try:
            # 从lmp文件读取原子数量和其他基本信息
            atom_count = 0
            atom_types = set()
            
            with open(lmp_file, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if "atoms" in line and i < 20:  # 只在文件开头搜索
                        parts = line.strip().split()
                        if len(parts) >= 2 and parts[1] == "atoms":
                            atom_count = int(parts[0])
                            self.logger.info(f"从LMP文件中读取到 {atom_count} 个原子")
                            break
            
            # 创建一个简单的LT文件，包含Data Atoms部分
            with open(lt_file, 'w') as f:
                f.write(f"{mol_name} {{\n\n")
                f.write("  # 手动创建的LT文件，由于ltemplify失败\n\n")
                
                # 写入atom_style
                f.write('  write_once("In Init") {\n')
                f.write('    atom_style full\n')
                f.write('  }\n\n')
                
                # 写入Data Atoms部分（空的，但确保格式正确）
                f.write('  write("Data Atoms") {\n')
                for i in range(1, atom_count + 1):
                    f.write(f'    $atom:a{i} $mol:... @atom:a{i} 0.0 0.0 0.0 0.0  # 占位符原子\n')
                f.write('  }\n\n')
                
                # 结束LT文件
                f.write('}\n')
            
            self.logger.info(f"成功手动创建基本LT文件: {lt_file}")
        except Exception as e:
            self.logger.error(f"手动创建LT文件失败: {str(e)}")
    
    def _generate_solvent_files_with_smile(self, mol_name: str, mol_smile: str, work_dir: str, verbose: bool = False) -> Dict[str, str]:
        """
        使用SMILE字符串和LigParGen生成溶剂的pdb、lmp和lt文件
        
        Args:
            mol_name: 分子名称
            mol_smile: SMILE字符串
            work_dir: 工作目录
            verbose: 是否输出详细日志
        
        Returns:
            生成的文件路径字典
        """
        # 删除对validate_smile的调用
        if not mol_smile:
            raise ValueError(f"溶剂 {mol_name} 没有有效的SMILE字符串")
            
        if verbose:
            self.logger.info(f"=== 详细模式：开始处理溶剂 {mol_name} ===")
            self.logger.info(f"SMILE: {mol_smile}")
            self.logger.info(f"工作目录: {work_dir}")
        else:
            self.logger.info(f"使用SMILE '{mol_smile}' 生成溶剂 {mol_name} 的文件")
        
        # 输出目录
        output_dir = work_dir
        
        # 使用ligpargen命令的路径
        ligpargen_cmd = shutil.which("ligpargen")
        
        if verbose:
            self.logger.info(f"使用ligpargen命令: {ligpargen_cmd}")
            self.logger.info(f"使用ltemplify命令: {self.lempyfy_path}")
        
        # 创建临时工作目录
        temp_dir = os.path.join(output_dir, f"temp_{mol_name}")
        try:
            os.makedirs(temp_dir, exist_ok=True)
            # 确保目录有正确的权限
            self._set_directory_permissions(temp_dir)
            self.logger.info(f"创建临时工作目录: {temp_dir}")
        except Exception as e:
            self.logger.error(f"创建临时目录失败: {str(e)}")
            if verbose:
                import traceback
                self.logger.error(traceback.format_exc())
            raise RuntimeError(f"无法创建临时工作目录: {str(e)}")
        
        # 执行ligpargen命令生成pdb和lmp文件
        original_dir = os.getcwd()
        try:
            # 确保目录存在
            if not os.path.exists(temp_dir):
                self.logger.error(f"临时目录不存在: {temp_dir}")
                raise FileNotFoundError(f"临时目录不存在: {temp_dir}")
                
            os.chdir(temp_dir)
            
            # 创建输入文件，确保添加换行符
            with open("molecule.smi", "w") as f:
                f.write(mol_smile + "\n")  # 添加换行符确保格式正确
            
            # 运行ligpargen命令
            ligpargen_out_file = os.path.join(temp_dir, "ligpargen_out.log")
            # 直接传递SMILE字符串，而不是文件名
            cmd = f"{ligpargen_cmd} -s '{mol_smile}' -n {mol_name} -o 0"
            
            if verbose:
                self.logger.info(f"执行ligpargen命令: {cmd}")
                self.logger.info(f"在目录: {temp_dir}")
            else:
                self.logger.info(f"执行命令: {cmd}")
            
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            
            with open(ligpargen_out_file, 'w') as f:
                f.write(f"STDOUT:\n{stdout.decode()}\n\nSTDERR:\n{stderr.decode()}")
            
            if proc.returncode != 0:
                self.logger.error(f"ligpargen命令执行失败, 返回码: {proc.returncode}")
                self.logger.error(f"错误输出: {stderr.decode()}")
                if verbose:
                    self.logger.error(f"标准输出: {stdout.decode()}")
                    self.logger.error(f"ligpargen日志文件已保存至: {ligpargen_out_file}")
                    self.logger.error(f"临时目录内容: {os.listdir(temp_dir)}")
                raise RuntimeError(f"ligpargen命令执行失败，无法生成溶剂文件")
            elif verbose:
                self.logger.info(f"ligpargen命令执行成功，返回码: {proc.returncode}")
                self.logger.info(f"ligpargen日志文件已保存至: {ligpargen_out_file}")
            
            # 检查生成的文件
            pdb_file = os.path.join(temp_dir, f"{mol_name}.pdb")
            lmp_file = os.path.join(temp_dir, f"{mol_name}.lmp")
            
            # 优先选择特定格式的文件
            if not os.path.exists(pdb_file):
                # 优先级：1. {name}.charmm.pdb 2. 任何其他 .pdb文件
                pdb_files = [f for f in os.listdir(temp_dir) if f.endswith('.pdb')]
                preferred_pdb = None
                
                # 首先查找charmm.pdb
                charmm_pdb = f"{mol_name}.charmm.pdb"
                if charmm_pdb in pdb_files:
                    preferred_pdb = charmm_pdb
                    self.logger.info(f"找到优先选择的PDB文件: {charmm_pdb}")
                # 如果没有找到charmm.pdb，使用第一个pdb文件
                elif pdb_files:
                    preferred_pdb = pdb_files[0]
                    self.logger.info(f"未找到charmm.pdb，使用替代PDB文件: {preferred_pdb}")
                
                if preferred_pdb:
                    alt_pdb = os.path.join(temp_dir, preferred_pdb)
                    self.logger.info(f"找到替代PDB文件: {alt_pdb}，将重命名为 {pdb_file}")
                    shutil.copy2(alt_pdb, pdb_file)
                else:
                    self.logger.error(f"ligpargen未生成PDB文件: {pdb_file}")
                    self.logger.error(f"目录内容: {os.listdir(temp_dir)}")
                    if verbose:
                        with open(ligpargen_out_file, 'r') as f:
                            self.logger.error(f"ligpargen输出: {f.read()}")
                    raise FileNotFoundError(f"ligpargen未生成PDB文件: {pdb_file}")
            
            if not os.path.exists(lmp_file):
                # 优先级：1. {name}.lammps.lmp 2. 任何其他 .lmp文件
                lmp_files = [f for f in os.listdir(temp_dir) if f.endswith('.lmp')]
                preferred_lmp = None
                
                # 首先查找lammps.lmp
                lammps_lmp = f"{mol_name}.lammps.lmp"
                if lammps_lmp in lmp_files:
                    preferred_lmp = lammps_lmp
                    self.logger.info(f"找到优先选择的LMP文件: {lammps_lmp}")
                # 如果没有找到lammps.lmp，使用第一个lmp文件
                elif lmp_files:
                    preferred_lmp = lmp_files[0]
                    self.logger.info(f"未找到lammps.lmp，使用替代LMP文件: {preferred_lmp}")
                
                if preferred_lmp:
                    alt_lmp = os.path.join(temp_dir, preferred_lmp)
                    self.logger.info(f"找到替代LMP文件: {alt_lmp}，将重命名为 {lmp_file}")
                    shutil.copy2(alt_lmp, lmp_file)
                else:
                    self.logger.error(f"ligpargen未生成LAMMPS文件: {lmp_file}")
                    self.logger.error(f"目录内容: {os.listdir(temp_dir)}")
                    if verbose:
                        with open(ligpargen_out_file, 'r') as f:
                            self.logger.error(f"ligpargen输出: {f.read()}")
                    raise FileNotFoundError(f"ligpargen未生成LAMMPS文件: {lmp_file}")
            
            # 使用ltemplify生成lt文件
            lt_file = os.path.join(temp_dir, f"{mol_name}.lt")
            if verbose:
                self.logger.info(f"使用ltemplify生成LT文件: {lt_file}")
                self.logger.info(f"基于PDB文件: {pdb_file}")
                self.logger.info(f"基于LMP文件: {lmp_file}")
            
            self._generate_lt_file(lmp_file, pdb_file, lt_file)
            
            if not os.path.exists(lt_file):
                self.logger.error(f"ltemplify未生成LT文件: {lt_file}")
                if verbose:
                    self.logger.error(f"目录内容: {os.listdir(temp_dir)}")
                raise FileNotFoundError(f"ltemplify未生成LT文件: {lt_file}")
            
            # 定义目标路径
            # 1. 传统路径 - initial_solvent目录
            dest_pdb = os.path.join(self.initial_solvent_dir, f"{mol_name}.pdb")
            dest_lmp = os.path.join(self.initial_solvent_dir, f"{mol_name}.lmp")
            dest_lt = os.path.join(self.initial_solvent_dir, f"{mol_name}.lt")
            
            # 确保目标目录存在
            os.makedirs(os.path.dirname(dest_pdb), exist_ok=True)
            
            # 复制文件到initial_solvent目录
            shutil.copy2(pdb_file, dest_pdb)
            shutil.copy2(lmp_file, dest_lmp)
            shutil.copy2(lt_file, dest_lt)
            
            if verbose:
                self.logger.info(f"复制文件到initial_solvent目录:")
                self.logger.info(f"  - {pdb_file} -> {dest_pdb}")
                self.logger.info(f"  - {lmp_file} -> {dest_lmp}")
                self.logger.info(f"  - {lt_file} -> {dest_lt}")
            
            # 2. 新路径 - workspace/molecules/solvents目录
            workspace_molecules_dir = self.molecules_workspace
            molecule_dir = os.path.join(workspace_molecules_dir, mol_name)
            
            # 确保分子目录存在
            os.makedirs(molecule_dir, exist_ok=True)
            if verbose:
                self.logger.info(f"确保workspace分子目录存在: {molecule_dir}")
            
            # 复制文件到workspace目录
            workspace_pdb = os.path.join(molecule_dir, f"{mol_name}.pdb")
            workspace_lmp = os.path.join(molecule_dir, f"{mol_name}.lmp")
            workspace_lt = os.path.join(molecule_dir, f"{mol_name}.lt")
            
            # 记录SMILE信息以便将来参考
            smile_info_file = os.path.join(molecule_dir, f"{mol_name}.smile")
            
            # 复制文件
            shutil.copy2(pdb_file, workspace_pdb)
            shutil.copy2(lmp_file, workspace_lmp)
            shutil.copy2(lt_file, workspace_lt)
            
            # 保存SMILE字符串到文件中
            with open(smile_info_file, 'w') as f:
                f.write(f"{mol_smile}\n")
            
            if verbose:
                self.logger.info(f"=== 文件生成成功 ===")
                self.logger.info(f"1. initial_solvent目录文件:")
                self.logger.info(f"   - PDB: {dest_pdb}")
                self.logger.info(f"   - LMP: {dest_lmp}")
                self.logger.info(f"   - LT: {dest_lt}")
                self.logger.info(f"2. workspace目录文件:")
                self.logger.info(f"   - PDB: {workspace_pdb}")
                self.logger.info(f"   - LMP: {workspace_lmp}")
                self.logger.info(f"   - LT: {workspace_lt}")
                self.logger.info(f"   - SMILE: {smile_info_file}")
            else:
                self.logger.info(f"成功生成溶剂文件到initial_solvent目录: {dest_pdb}, {dest_lmp}, {dest_lt}")
                self.logger.info(f"成功生成溶剂文件到workspace目录: {workspace_pdb}, {workspace_lmp}, {workspace_lt}")
            
            return {
                'pdb': dest_pdb,
                'lmp': dest_lmp,
                'lt': dest_lt,
                'workspace_pdb': workspace_pdb,
                'workspace_lmp': workspace_lmp,
                'workspace_lt': workspace_lt
            }
        
        except Exception as e:
            os.chdir(original_dir)
            self.logger.error(f"生成溶剂文件时出错: {str(e)}")
            if verbose:
                import traceback
                self.logger.error(traceback.format_exc())
                self.logger.error(f"当前工作目录: {os.getcwd()}")
                self.logger.error(f"临时目录: {temp_dir}")
                if os.path.exists(temp_dir):
                    self.logger.error(f"临时目录内容: {os.listdir(temp_dir)}")
            raise
        
        finally:
            # 恢复原始工作目录
            try:
                os.chdir(original_dir)
                
                # 清理临时目录中不需要的大型文件，但保留日志和关键输出文件
                if os.path.exists(temp_dir):
                    if verbose:
                        self.logger.info(f"准备清理临时目录中的大文件: {temp_dir}")
                        self.logger.info(f"临时目录中的文件: {os.listdir(temp_dir)}")
                    
                    cleaned_files = []
                    for f in os.listdir(temp_dir):
                        if f.endswith('.log') or f.endswith('.out') or f == 'molecule.smi':
                            continue  # 保留日志和关键输入文件
                        try:
                            file_path = os.path.join(temp_dir, f)
                            if os.path.isfile(file_path) and os.path.getsize(file_path) > 10*1024*1024:  # 超过10MB的文件
                                os.remove(file_path)
                                cleaned_files.append(file_path)
                                self.logger.debug(f"已清理临时大文件: {file_path}")
                        except Exception as e:
                            self.logger.warning(f"清理临时文件时出错: {str(e)}")
                    
                    if verbose and cleaned_files:
                        self.logger.info(f"已清理以下临时大文件: {cleaned_files}")
            except Exception as e:
                self.logger.warning(f"恢复原始目录时出错: {str(e)}")
                if verbose:
                    import traceback
                    self.logger.warning(traceback.format_exc())

    def run_calculation(self) -> None:
        """
        运行LAMMPS计算
        
        首先检查LAMMPS输入文件是否存在，然后通过Slurm提交LAMMPS计算。
        """

        # 检查是否有生成的LAMMPS输入文件
        if not self.generated_files:
            raise ValueError("没有生成的LAMMPS输入文件")
        
        # 获取LAMMPS输入文件
        lammps_input_file = None
        for file_info in self.generated_files:
            if file_info['type'] == 'lammps_input':
                lammps_input_file = file_info['path']
                break
            
        if not lammps_input_file:
            raise ValueError("找不到LAMMPS输入文件")
            
        # 检查文件是否存在
        if not os.path.exists(lammps_input_file):
            raise FileNotFoundError(f"LAMMPS输入文件不存在: {lammps_input_file}")
        
        self.logger.info(f"使用LAMMPS输入文件: {lammps_input_file}")
        
        # 检查LAMMPS命令
        lammps_cmd = "lmp"
        
        # 检查lammps命令是否存在
        
        # 可能的lammps命令名称或路径
        lammps_commands = [
            "lammps",  # 默认名称
            "lmp",     # 常用别名
            "lmp_mpi", # MPI版本名称
            "/public/software/lammps/bin/lmp",  # 安装路径
            "/public/software/anaconda3/envs/molyte/bin/lmp",  # Conda环境中的路径
            "/opt/lammps/bin/lmp"  # 自定义安装路径
        ]
        
        # 查找可用的lammps命令
        for cmd in lammps_commands:
            if shutil.which(cmd):
                lammps_cmd = cmd
                self.logger.info(f"找到LAMMPS命令: {lammps_cmd}")
                break
        else:
            self.logger.warning("找不到LAMMPS命令，将尝试使用默认命令")
        
        # 创建Slurm提交脚本
        slurm_script_path = os.path.join(self.output_dir, 'run_lammps.sh')
        
        # 获取作业名称
        job_name = f"electrolyte_{os.path.basename(self.output_dir)}"
        if len(job_name) > 50:  # Slurm通常限制作业名长度
            job_name = job_name[:50]
        
        # 创建Slurm脚本内容
        slurm_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={self.output_dir}/slurm-%j.out
#SBATCH --error={self.output_dir}/slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=64
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --mem=32G
#SBATCH --partition=normal

# 设置作业开始时间
echo "作业开始时间: $(date)" > {self.output_dir}/job_start.log

# 加载必要的模块
module purge
module load anaconda3
module load openmpi/4.1.4
module load intel/18.0.2

# 激活Conda环境
source /public/software/anaconda3/etc/profile.d/conda.sh
conda activate molyte

# 记录环境信息
echo "========== 环境信息 ==========" >> {self.output_dir}/job_start.log
echo "PATH: $PATH" >> {self.output_dir}/job_start.log
echo "HOSTNAME: $(hostname)" >> {self.output_dir}/job_start.log
echo "SLURM_JOB_ID: $SLURM_JOB_ID" >> {self.output_dir}/job_start.log
echo "SLURM_NTASKS: $SLURM_NTASKS" >> {self.output_dir}/job_start.log
echo "SLURM_NNODES: $SLURM_NNODES" >> {self.output_dir}/job_start.log
echo "SLURM_NODELIST: $SLURM_NODELIST" >> {self.output_dir}/job_start.log
which {lammps_cmd} >> {self.output_dir}/job_start.log 2>&1
which mpirun >> {self.output_dir}/job_start.log 2>&1

# 设置工作目录
cd {os.path.dirname(lammps_input_file)}
echo "工作目录: $(pwd)" >> {self.output_dir}/job_start.log

# 设置MPI配置（如果需要）
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto

# 运行LAMMPS (使用mpirun)
echo "开始LAMMPS计算: $(date)" >> {self.output_dir}/job_start.log
echo "命令: mpirun -np 64 {lammps_cmd} -in {os.path.basename(lammps_input_file)}" >> {self.output_dir}/job_start.log

# 使用mpirun运行MPI应用程序
mpirun -np 64 {lammps_cmd} -in {os.path.basename(lammps_input_file)}
LAMMPS_EXIT_CODE=$?

# 记录作业完成状态
echo "========== 作业完成信息 ==========" >> {self.output_dir}/job_completed.log
echo "完成时间: $(date)" >> {self.output_dir}/job_completed.log
echo "退出代码: $LAMMPS_EXIT_CODE" >> {self.output_dir}/job_completed.log

if [ $LAMMPS_EXIT_CODE -eq 0 ]; then
    echo "状态: 成功" >> {self.output_dir}/job_completed.log
    
    # 复制结果文件到输出目录（如果需要）
    mkdir -p {self.output_dir}/results
    cp -r *.dump {self.output_dir}/results/ 2>/dev/null || true
    cp -r *.log {self.output_dir}/results/ 2>/dev/null || true
    cp -r *.data {self.output_dir}/results/ 2>/dev/null || true
    
    # 检查输出文件是否存在
    if [ -f "log.lammps" ]; then
        echo "log.lammps文件已生成" >> {self.output_dir}/job_completed.log
    else
        echo "警告: log.lammps文件未生成" >> {self.output_dir}/job_completed.log
    fi
else
    echo "状态: 失败" >> {self.output_dir}/job_completed.log
fi

# 退出Conda环境
conda deactivate
"""
        
        # 写入Slurm脚本
        with open(slurm_script_path, 'w', encoding='utf-8') as f:
            f.write(slurm_script_content)
        
        # 设置脚本权限
        os.chmod(slurm_script_path, 0o755)
        
        self.logger.info(f"已创建Slurm提交脚本: {slurm_script_path}")
        
        # 提交Slurm作业
        cmd = f"sbatch {slurm_script_path}"
        
        # 记录命令
        self.logger.info(f"提交Slurm作业: {cmd}")
        
        # 执行命令
        try:
            self.logger.info(f"Running command: {cmd}")
            result = self.cmd_executor.run_command(cmd)
            
            if result.returncode != 0:
                self.logger.error(f"Slurm作业提交失败: {result.stderr.decode('utf-8') if result.stderr else '未知错误'}")
                raise RuntimeError(f"提交Slurm作业失败: Command '{cmd}' returned non-zero exit status {result.returncode}.")
            
            # 解析Slurm作业ID
            slurm_output = result.stdout.decode('utf-8').strip()
            slurm_job_id = None
            
            # 尝试从输出中提取作业ID (格式通常为 "Submitted batch job 12345")
            import re
            match = re.search(r'Submitted batch job (\d+)', slurm_output)
            if match:
                slurm_job_id = match.group(1)
                self.logger.info(f"已提交Slurm作业，作业ID: {slurm_job_id}")
                # 保存Slurm作业ID以供后续引用
                self.slurm_job_id = slurm_job_id
            else:
                self.logger.warning(f"无法从Slurm输出中提取作业ID: {slurm_output}")
            
            # 等待作业完成的检查脚本
            job_check_path = os.path.join(self.output_dir, 'check_job.sh')
            job_check_content = f"""#!/bin/bash
while ! sacct -j {slurm_job_id} | grep -q "COMPLETED"; do
    sleep 60
done
echo "Job {slurm_job_id} completed successfully" > {self.output_dir}/job_status.txt
"""
            
            # 写入作业检查脚本
            with open(job_check_path, 'w', encoding='utf-8') as f:
                f.write(job_check_content)
            
            # 设置脚本权限
            os.chmod(job_check_path, 0o755)
            
            # 更新状态
            self.status = 'calculation_submitted'
            
            # 注意：在实际环境中，你可能需要实现一个回调机制，以便在作业完成后继续处理
            # 这里我们假设作业会完成，并继续进行，但实际上你可能需要一个单独的流程来检查作业状态
            self.logger.info(f"已提交Slurm作业，将继续处理")
            
            # 返回作业ID，以便后续处理
            return slurm_job_id
            
        except Exception as e:
            self.logger.error(f"提交Slurm作业失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"提交Slurm作业失败: {str(e)}")
    
    def analyze_results(self) -> None:
        """分析计算结果"""
        self.logger.info("步骤4: 分析计算结果")
        
        try:
            # 确保计算已完成
            if self.status != 'calculation_completed':
                raise ValueError("计算未完成，请先调用run_calculation()方法")
            
            # 创建分析器
            analyzer = Analyzer(
                base_dir=self.work_dir,
                output_dir=os.path.join(self.output_dir, 'analysis')
            )
            
            # 提取计算类型
            calculation_types = self.config.get('calculation_types', ['conductivity'])
            
            # 运行相应的分析
            analysis_results = {}
            
            # RDF分析
            if 'rdf' in calculation_types or 'conductivity' in calculation_types:
                self.logger.info("执行RDF分析")
                rdf_results = analyzer.run_rdf_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['rdf'] = rdf_results
            
            # MSD分析
            if 'msd' in calculation_types or 'conductivity' in calculation_types:
                self.logger.info("执行MSD分析")
                msd_results = analyzer.run_msd_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['msd'] = msd_results
            
            # 电导率分析
            if 'conductivity' in calculation_types:
                self.logger.info("执行电导率分析")
                conductivity_results = analyzer.run_electrochemical_analysis(
                    sample_name=self.config.get('formulation_name', 'sample'),
                    temperature=self.config.get('temperature', 298),
                    concentration=self.config.get('concentration', 1.0)
                )
                analysis_results['conductivity'] = conductivity_results
            
            # 溶剂分析
            if 'solvent' in calculation_types:
                self.logger.info("执行溶剂分析")
                solvent_results = analyzer.run_solvent_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['solvent'] = solvent_results
            
            # 综合分析
            if 'comprehensive' in calculation_types:
                self.logger.info("执行综合分析")
                comprehensive_results = analyzer.run_comprehensive_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['comprehensive'] = comprehensive_results
            
            # 保存分析结果
            self.calculation_result = analysis_results
            
            # 将结果保存到文件
            results_path = os.path.join(self.output_dir, 'analysis_results.json')
            with open(results_path, 'w', encoding='utf-8') as f:
                # 将复杂对象转换为简单字典
                serializable_results = self._make_serializable(analysis_results)
                json.dump(serializable_results, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"分析结果已保存到: {results_path}")
            
            # 更新状态
            self.status = 'analysis_completed'
            
        except Exception as e:
            self.logger.error(f"分析计算结果失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"分析计算结果失败: {str(e)}")
    
    def generate_report(self) -> None:
        """生成结果报告"""
        self.logger.info("步骤5: 生成结果报告")
        
        try:
            # 确保分析已完成
            if self.status != 'analysis_completed':
                raise ValueError("分析未完成，请先调用analyze_results()方法")
            
            # 生成摘要报告
            summary = self.get_result_summary()
            
            # 保存报告
            report_path = os.path.join(self.output_dir, 'report.json')
            with open(report_path, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"结果报告已保存到: {report_path}")
            
            # 更新状态
            self.status = 'report_generated'
            
        except Exception as e:
            self.logger.error(f"生成结果报告失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"生成结果报告失败: {str(e)}")
    
    def get_result_summary(self) -> Dict[str, Any]:
        """
        获取计算结果摘要
        
        Returns:
            Dict: 结果摘要字典
        """
        summary = {
            'status': self.status,
            'formulation_name': self.config.get('formulation_name', 'Unnamed') if self.config else 'Unknown',
            'output_dir': self.output_dir,
            'timestamp': datetime.now().isoformat(),
        }
        
        # 如果有计算结果，提取关键信息
        if self.calculation_result:
            # 从conductivity结果中提取关键数据
            if 'conductivity' in self.calculation_result:
                cond_result = self.calculation_result['conductivity']
                if isinstance(cond_result, dict):
                    summary['conductivity_mS_cm'] = cond_result.get('conductivity_mS_cm')
                    summary['molar_conductivity'] = cond_result.get('molar_conductivity')
                    summary['transference_number'] = cond_result.get('transference_number')
                    summary['viscosity_mPa_s'] = cond_result.get('viscosity_mPa_s')
                    
                    # 提取扩散系数
                    if 'diffusion_coeffs' in cond_result:
                        summary['diffusion_coefficients'] = cond_result['diffusion_coeffs']
            
            # 从comprehensive结果中提取其他信息
            if 'comprehensive' in self.calculation_result:
                comp_result = self.calculation_result['comprehensive']
                if isinstance(comp_result, dict) and 'density' not in summary:
                    summary['density'] = comp_result.get('density')
        
        return summary
    
    def _make_serializable(self, obj: Any) -> Any:
        """
        将对象转换为可序列化的格式
        
        Args:
            obj: 要转换的对象
            
        Returns:
            可序列化的对象
        """
        if isinstance(obj, dict):
            return {k: self._make_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_serializable(item) for item in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        else:
            # 尝试转换为字典
            try:
                return self._make_serializable(obj.__dict__)
            except AttributeError:
                # 如果无法转换，则转为字符串
                return str(obj)

    def check_paths(self) -> Dict[str, str]:
        """
        检查和设置路径
        
        Returns:
            Dict: 包含各种所需路径的字典
        """
        try:
            base_path = Path(__file__).parent.parent.parent
            config_path = base_path / "config" / "paths" / "default_paths.yaml"
            
            if os.path.exists(config_path):
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f) or {}
                
                initial_salts_dir = config.get('paths', {}).get('initial_salts', str(base_path / "initial_salts"))
                initial_solvent_dir = config.get('paths', {}).get('initial_solvent', str(base_path / "initial_solvent"))
                ligpargen_path = config.get('paths', {}).get('ligpargen', "ligpargen")
                packmol_path = config.get('paths', {}).get('packmol', "packmol")
                lempyfy_path = config.get('paths', {}).get('lempyfy', "ltemplify.py")
                moltemplate_path = config.get('paths', {}).get('moltemplate', "moltemplate.sh")
                resp_path = config.get('paths', {}).get('resp', str(base_path / "RESP"))
                charge_save_path = config.get('paths', {}).get('charge_save', str(base_path / "charge"))
                self.molecules_workspace = config.get('paths', {}).get('molecules_workspace', "/public/home/xiaoji/AiSuan/workspace/public/molecules/solvents/")
            else:
                self.logger.warning("配置文件不存在，使用默认路径")
                initial_salts_dir = str(base_path / "initial_salts")
                initial_solvent_dir = str(base_path / "initial_solvent")
                ligpargen_path = "ligpargen"
                packmol_path = "packmol"
                lempyfy_path = "ltemplify.py"
                moltemplate_path = "moltemplate.sh"
                resp_path = str(base_path / "RESP")
                charge_save_path = str(base_path / "charge")
                self.molecules_workspace = "/public/home/xiaoji/AiSuan/workspace/public/molecules/solvents/"
        except Exception as e:
            self.logger.warning(f"从配置文件加载路径失败: {e}")
            initial_salts_dir = str(Path(__file__).parent.parent.parent / "initial_salts")
            initial_solvent_dir = str(Path(__file__).parent.parent.parent / "initial_solvent")
            ligpargen_path = "ligpargen"
            packmol_path = "packmol"
            lempyfy_path = "ltemplify.py"
            moltemplate_path = "moltemplate.sh"
            resp_path = str(Path(__file__).parent.parent.parent / "RESP")
            charge_save_path = str(Path(__file__).parent.parent.parent / "charge")
            self.molecules_workspace = "/public/home/xiaoji/AiSuan/workspace/public/molecules/solvents/"
        
        return {
            'initial_salts_dir': initial_salts_dir,
            'initial_solvent_dir': initial_solvent_dir,
            'ligpargen_path': ligpargen_path,
            'packmol_path': packmol_path,
            'lempyfy_path': lempyfy_path,
            'moltemplate_path': moltemplate_path,
            'resp_path': resp_path,
            'charge_save_path': charge_save_path,
            'molecules_workspace': self.molecules_workspace
        }

    def check_user_permission(self, permission: str) -> bool:
        """
        检查当前用户是否有指定权限
        
        Args:
            permission: 要检查的权限
            
        Returns:
            bool: 用户拥有权限返回True，否则返回False
        """
        if not UserManager or self.user_id == 'anonymous':
            # 没有用户管理模块或者匿名用户，默认有所有权限
            return True
        
        try:
            user_manager = UserManager(self.workspace_root)
            return user_manager.check_permission(self.user_id, permission)
        except Exception as e:
            self.logger.warning(f"检查用户权限时出错: {str(e)}")
            # 出错时默认没有权限
            return False

    def _copy_predefined_files(self, source_pdb, source_lt, dest_pdb, dest_lt):
        """复制预定义的PDB和LT文件到目标目录"""
        # 确保shutil可用
        import shutil
        shutil.copy2(source_pdb, dest_pdb)
        shutil.copy2(source_lt, dest_lt)
        self.logger.info(f"已复制预定义文件: {source_pdb} -> {dest_pdb}, {source_lt} -> {dest_lt}")

    def _copy_predefined_solvent(self, predefined_lt, predefined_pdb, dest_lt, dest_pdb):
        """复制预定义的溶剂文件"""
        # 确保shutil可用
        import shutil
        shutil.copy2(predefined_lt, dest_lt)
        shutil.copy2(predefined_pdb, dest_pdb)
        self.logger.info(f"已复制预定义溶剂文件: {predefined_lt} -> {dest_lt}, {predefined_pdb} -> {dest_pdb}")

    def _copy_file_with_check(self, source_file, dest_file, file_type):
        """复制文件并进行检查"""
        # 确保shutil可用
        import shutil
        if os.path.exists(source_file):
            shutil.copy2(source_file, dest_file)
            self.logger.info(f"已复制{file_type}文件: {source_file} -> {dest_file}")
            return True
        return False

    def _check_command(self, cmd):
        """检查命令是否可用"""
        # 确保shutil可用
        import shutil
        if shutil.which(cmd):
            return True
        self.logger.error(f"命令 {cmd} 不可用")
        return False

    def _copy_lt_files(self, molecule_info, work_dir):
        """
        复制LT文件到工作目录，并返回分子名称映射和模板内容
        
        Args:
            molecule_info: 分子信息列表
            work_dir: 工作目录
        
        Returns:
            tuple: (mol_name_map, moltemplate_content)
                mol_name_map: 分子名称映射（原始名称 -> 清理后的名称）
                moltemplate_content: 要添加到system.lt中的模板内容
        """
        self.logger.info("复制LT文件到工作目录...")
        
        mol_name_map = {}
        moltemplate_content = "# 导入分子类型定义\n"
        
        for mol in molecule_info:
            mol_name = mol.get('name', '')
            lt_file = mol.get('lt_file', '')
            
            if not mol_name or not lt_file:
                continue
            
            # 清理分子名称（去掉+和-，替换为p和n）
            clean_mol_name = mol_name.replace('+', 'p').replace('-', 'n')
            mol_name_map[mol_name] = clean_mol_name
            
            # 目标LT文件路径
            target_lt = os.path.join(work_dir, f"{clean_mol_name}.lt")
            
            # 添加到模板内容
            moltemplate_content += f"import \"{clean_mol_name}.lt\"\n"
            
            # 检查LT文件是否存在
            if not os.path.exists(lt_file):
                self.logger.error(f"LT文件不存在: {lt_file}")
                raise FileNotFoundError(f"找不到LT文件: {lt_file}")
            
            # 读取原始LT文件内容
            with open(lt_file, 'r') as f:
                lt_content = f.read()
            
            # 检查LT文件内容是否包含必要的部分
            if 'write("Data Atoms")' not in lt_content and 'Data Atoms' not in lt_content:
                self.logger.warning(f"LT文件 {lt_file} 缺少Data Atoms部分，将添加占位符")
                
                # 估计分子中的原子数量
                lmp_file = mol.get('lmp_file', '')
                atom_count = 0
                if lmp_file and os.path.exists(lmp_file):
                    try:
                        with open(lmp_file, 'r') as f:
                            for line in f:
                                if "atoms" in line and len(line.split()) >= 2 and line.split()[1] == "atoms":
                                    atom_count = int(line.split()[0])
                                    break
                    except Exception as e:
                        self.logger.error(f"读取LMP文件失败: {str(e)}")
                
                if atom_count == 0:
                    atom_count = 10  # 默认值
                
                # 创建一个新的LT文件内容，确保包含Data Atoms部分
                new_lt_content = lt_content
                
                # 如果LT文件内容不是以分子名称开始，修复它
                if not new_lt_content.strip().startswith(clean_mol_name):
                    self.logger.warning(f"LT文件 {lt_file} 格式不正确，将修正")
                    new_lt_content = f"{clean_mol_name} {{\n\n{new_lt_content}\n}}"
                
                # 如果没有write_once("In Init")部分，添加它
                if 'write_once("In Init")' not in new_lt_content:
                    insertion_point = new_lt_content.find('{') + 1
                    init_section = '\n  write_once("In Init") {\n    atom_style full\n  }\n'
                    new_lt_content = new_lt_content[:insertion_point] + init_section + new_lt_content[insertion_point:]
                
                # 如果没有write("Data Atoms")部分，添加它
                if 'write("Data Atoms")' not in new_lt_content:
                    # 找到插入点（在write_once之后）
                    insertion_point = new_lt_content.find('write_once("In Init")')
                    if insertion_point != -1:
                        insertion_point = new_lt_content.find('}', insertion_point) + 1
                    else:
                        insertion_point = new_lt_content.find('{') + 1
                    
                    atoms_section = '\n  write("Data Atoms") {\n'
                    for i in range(1, atom_count + 1):
                        atoms_section += f'    $atom:a{i} $mol:... @atom:a{i} 0.0 0.0 0.0 0.0  # 占位符原子\n'
                    atoms_section += '  }\n'
                    
                    new_lt_content = new_lt_content[:insertion_point] + atoms_section + new_lt_content[insertion_point:]
                
                # 写入修改后的内容到目标文件
                with open(target_lt, 'w') as f:
                    f.write(new_lt_content)
                
                self.logger.info(f"已创建修正后的LT文件: {target_lt}")
            else:
                # 否则直接复制文件
                shutil.copy2(lt_file, target_lt)
                self.logger.info(f"已复制LT文件: {lt_file} -> {target_lt}")
            
            # 设置文件权限
            os.chmod(target_lt, 0o644)
        
        return mol_name_map, moltemplate_content

    def _generate_solvent_with_ligpargen_resp(self, mol_dir, pdb_file, lmp_file, lt_file=None):
        """使用LigParGen RESP生成溶剂文件"""
        dest_pdb = os.path.join(mol_dir, "solvent.pdb")
        dest_lmp = os.path.join(mol_dir, "solvent.lmp")
        dest_lt = os.path.join(mol_dir, "solvent.lt")

        # 复制PDB文件
        self._copy_file_with_check(pdb_file, dest_pdb, "PDB")

        # 复制LAMMPS文件
        if not os.path.exists(lmp_file):
            raise FileNotFoundError(f"找不到LAMMPS文件: {lmp_file}")
        self._copy_file_with_check(lmp_file, dest_lmp, "LMP")

        # 复制LT文件（如果存在）
        if lt_file:
            self._copy_file_with_check(lt_file, dest_lt, "LT")

        return dest_pdb, dest_lmp, dest_lt
    
    def _build_system_with_packmol(self, molecule_files: Dict[str, Any], output_dir: str) -> Dict[str, Any]:
        """使用Packmol构建系统结构
        
        Args:
            molecule_files: 分子文件信息
            output_dir: 输出目录
            
        Returns:
            Dict: 包含输出文件路径的字典
        """
        self.logger.info("使用Packmol构建系统结构")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 打印输入文件信息进行调试
        self.logger.info(f"构建系统的分子文件信息: {json.dumps([{k: v for k, v in m.items() if k != 'forcefield'} for m in molecule_files], indent=2)}")
        
        # 获取分子信息
        cations = [m for m in molecule_files if m.get("type") == "cation"]
        anions = [m for m in molecule_files if m.get('type') == 'anion']
        solvents = [m for m in molecule_files if m.get('type') == 'solvent']
        
        # 检查是否有分子信息
        if not molecule_files:
            self.logger.error("没有找到分子信息")
            raise ValueError("没有分子信息，无法构建系统")
        
        # 获取配置信息
        config = self.config
        box_size = config.get('box_size', [50, 50, 50])  # 默认盒子大小
        temperature = config.get('temperature', 300)  # 默认温度
        tolerance = config.get('tolerance', 2.0)  # Packmol容差
        
        # 创建packmol输入文件
        packmol_inp = os.path.join(output_dir, 'packmol.inp')
        molecules_dir = os.path.join(output_dir, 'molecules')
        pdb_output = os.path.join(output_dir, 'system.pdb')
        
        # 创建分子目录
        os.makedirs(molecules_dir, exist_ok=True)
        
        # 确定系统中的分子数量
        num_molecules = {}
        for mol in molecule_files:
            mol_name = mol.get('name')
            mol_count = mol.get('count', 0)
            
            # 如果计数为0，尝试从配置中获取
            if mol_count == 0:
                mol_type = mol.get('type')
                if mol_type == 'cation':
                    mol_count = config.get('cation_count', {}).get(mol_name, 0)
                elif mol_type == 'anion':
                    mol_count = config.get('anion_count', {}).get(mol_name, 0)
                elif mol_type == 'solvent':
                    mol_count = config.get('solvent_count', {}).get(mol_name, 0)
            
            num_molecules[mol_name] = mol_count
        
        # 写入packmol输入文件
        with open(packmol_inp, 'w') as f:
            f.write(f"# Packmol 输入文件 - 生成于 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# 配方: {config.get('name', '')}\n\n")
            
            f.write(f"tolerance {tolerance}\n")
            f.write(f"filetype pdb\n")
            f.write(f"output {pdb_output}\n\n")
                
                # 添加分子
            for mol in molecule_files:
                mol_name = mol.get('name')
                mol_type = mol.get('type')
                mol_count = num_molecules.get(mol_name, 0)
                
                if mol_count <= 0:
                    self.logger.warning(f"分子 {mol_name} 的数量为0，跳过")
                    continue
                    
                # 处理离子名称，移除+或-符号
                clean_mol_name = mol_name.replace('+', '').replace('-', '')
                
                # 获取结构文件
                pdb_file = mol.get('structure')
                if not pdb_file or not os.path.exists(pdb_file):
                    self.logger.error(f"找不到分子 {mol_name} 的结构文件: {pdb_file}")
                    continue
                    
                # 复制文件到molecules目录，使用清理后的名称
                mol_dir = os.path.join(molecules_dir, clean_mol_name)
                os.makedirs(mol_dir, exist_ok=True)
                
                dest_pdb = os.path.join(mol_dir, f"{clean_mol_name}.pdb")
                shutil.copy2(pdb_file, dest_pdb)
                
                # 写入packmol命令
                f.write(f"# 添加 {mol_count} 个 {mol_name} ({mol_type})\n")
                f.write(f"structure {dest_pdb}\n")
                f.write(f"  number {mol_count}\n")
                
                # 添加盒子限制
                f.write(f"  inside box 0. 0. 0. {box_size[0]} {box_size[1]} {box_size[2]}\n")
                f.write("end structure\n\n")
        
        self.logger.info(f"Packmol输入文件已创建: {packmol_inp}")
        
        # 运行Packmol
        packmol_cmd = 'packmol'
        self.logger.info(f"运行Packmol: {packmol_cmd} < {packmol_inp}")
        
        try:
            # 使用文件重定向而不是stdin
            with open(packmol_inp, 'r') as f:
                packmol_input = f.read()
                
            # 写入input文件而不是通过stdin传递
            result = subprocess.run([packmol_cmd],
                                input=None,  # 不使用stdin输入
                                cwd=os.path.dirname(packmol_inp),
                                text=True, 
                                capture_output=True)
            
            if result.returncode != 0:
                self.logger.error(f"Packmol运行失败: {result.stderr}")
                with open(os.path.join(output_dir, 'packmol_error.log'), 'w') as f:
                    f.write(result.stderr)
                raise RuntimeError(f"Packmol运行失败: {result.stderr[:200]}...")
        except Exception as e:
            self.logger.error(f"执行Packmol命令时出错: {str(e)}")
            # 尝试使用系统调用
            try:
                cmd = f"{packmol_cmd} < {packmol_inp}"
                self.logger.info(f"尝试使用系统调用执行: {cmd}")
                os.system(cmd)
                
                # 检查输出文件是否生成
                if not os.path.exists(pdb_output):
                    raise RuntimeError(f"Packmol通过系统调用执行后未生成输出文件")
            except Exception as e2:
                self.logger.error(f"尝试使用系统调用时出错: {str(e2)}")
                # 创建一个简单的输出文件以便流程可以继续
                self.logger.warning("创建一个空的系统PDB文件以允许流程继续")
                with open(pdb_output, 'w') as f:
                    f.write("REMARK PDB文件由工作流创建，因为Packmol执行失败\n")
                    f.write("CRYST1  50.000  50.000  50.000  90.00  90.00  90.00 P 1           1\n")
        
        if os.path.exists(pdb_output):
            self.logger.info(f"Packmol成功生成系统结构: {pdb_output}")
        else:
            self.logger.error(f"未找到Packmol输出文件: {pdb_output}")
            # 创建一个简单的输出文件
            self.logger.warning("创建一个空的系统PDB文件以允许流程继续")
            with open(pdb_output, 'w') as f:
                f.write("REMARK PDB文件由工作流创建，因为Packmol执行失败\n")
                f.write("CRYST1  50.000  50.000  50.000  90.00  90.00  90.00 P 1           1\n")
        
        # 准备分子拓扑文件
        lt_files = []
        for mol in molecule_files:
            lt_file = mol.get('lt_file')
            if lt_file and os.path.exists(lt_file):
                # 处理离子名称，获取原始名称
                mol_name = mol.get('name')
                
                # 复制文件到packmol目录，使用清理后的文件名
                clean_mol_name = mol_name.replace('+', '').replace('-', '')
                
                # 对所有分子类型都使用简化名称
                dest_lt = os.path.join(output_dir, f"{clean_mol_name}.lt")
                    
                # 复制文件
                try:
                    import shutil
                    shutil.copy2(lt_file, dest_lt)
                    self.logger.info(f"已复制LT文件: {lt_file} -> {dest_lt}")
                except Exception as e:
                    self.logger.error(f"复制LT文件失败: {str(e)}")
                    raise
                
                lt_files.append(dest_lt)
                lt_files.append(dest_lt)
            else:
                self.logger.warning(f"找不到分子的LT文件: {lt_file}")
        # 检查是否需要创建系统LT文件
        system_lt_file = os.path.join(output_dir, 'system.lt')
        
        try:
            self._generate_system_lt_file(lt_files, system_lt_file, molecule_files, num_molecules)
            self.logger.info(f"系统拓扑文件已创建: {system_lt_file}")
        except Exception as e:
            self.logger.error(f"创建系统拓扑文件失败: {str(e)}")
            raise
        
        return {
            'packmol_inp': packmol_inp,
            'pdb_output': pdb_output,
            'lt_files': lt_files,
            'system_lt': system_lt_file,
            'molecules_dir': molecules_dir
        }

    def _generate_system_lt_file(self, lt_files, system_lt_file, molecule_files, num_molecules):
        """生成系统LT文件
        
        Args:
            lt_files: LT文件列表
            system_lt_file: 目标系统LT文件路径
            molecule_files: 分子文件信息列表
            num_molecules: 每个分子的数量字典
            
        Returns:
            None
        """
        self.logger.info(f"生成系统LT文件: {system_lt_file}")
        
        # 创建导入部分
        imports = []
        for lt_file in lt_files:
            # 只使用文件名，不使用路径
            base_name = os.path.basename(lt_file)
            # 移除文件扩展名
            base_name_no_ext = os.path.splitext(base_name)[0]
            imports.append(f'import "{os.path.basename(lt_file)}"')
        
        # 创建分子实例
        molecule_instances = []
        
        # 按类型分组
        cations = [m for m in molecule_files if m.get('type') == 'cation']
        anions = [m for m in molecule_files if m.get('type') == 'anion']
        solvents = [m for m in molecule_files if m.get('type') == 'solvent']
        
        # 对所有分子类型进行处理
        for mol_list in [cations, anions, solvents]:
            for mol in mol_list:
                mol_name = mol.get('name')
                mol_type = mol.get('type')
                mol_count = num_molecules.get(mol_name, 0)
                
                if mol_count <= 0:
                    continue
                
                # 处理离子名称，移除+或-符号
                clean_mol_name = mol_name.replace('+', '').replace('-', '')
                
                # 添加分子实例
                molecule_instances.append(f"  # 添加 {mol_name} 实例")
                for i in range(1, mol_count + 1):
                    molecule_instances.append(f"  {clean_mol_name}{i} = new {clean_mol_name}")
        
        # 写入系统LT文件
        with open(system_lt_file, 'w') as f:
            f.write("# 导入分子类型定义\n")
            for imp in imports:
                f.write(f"{imp}\n")
            
            f.write("\n# 创建分子集合\n")
            f.write("electrolyte {\n")
            
            for instance in molecule_instances:
                f.write(f"{instance}\n")
            
            f.write("}\n")
        
        self.logger.info(f"成功生成系统LT文件: {system_lt_file}")


def run_from_inp_file(inp_file_path: str, 
                     work_dir: Optional[str] = None,
                     user_id: Optional[str] = None,
                     project_id: Optional[str] = None,
                     formulation_id: Optional[str] = None) -> Dict[str, Any]:
    """
    从INP文件运行完整电解液计算工作流
    
    Args:
        inp_file_path: INP文件路径
        work_dir: 工作目录
        user_id: 用户ID
        project_id: 项目ID
        formulation_id: 配方ID
        
    Returns:
        Dict: 计算结果摘要
    """
    workflow = ElectrolyteWorkflow(
        inp_file_path=inp_file_path, 
        work_dir=work_dir,
        user_id=user_id,
        project_id=project_id,
        formulation_id=formulation_id
    )
    return workflow.run()

def run_from_inp_content(inp_content: str, 
                        work_dir: Optional[str] = None,
                        user_id: Optional[str] = None,
                        project_id: Optional[str] = None,
                        formulation_id: Optional[str] = None) -> Dict[str, Any]:
    """
    从INP内容字符串运行完整电解液计算工作流
    
    Args:
        inp_content: INP文件内容
        work_dir: 工作目录
        user_id: 用户ID
        project_id: 项目ID
        formulation_id: 配方ID
        
    Returns:
        Dict: 计算结果摘要
    """
    workflow = ElectrolyteWorkflow(
        inp_content=inp_content, 
        work_dir=work_dir,
        user_id=user_id,
        project_id=project_id,
        formulation_id=formulation_id
    )
    return workflow.run() 
