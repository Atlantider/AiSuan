"""
电解液工作流模块 - 整合版
该模块提供了电解液计算的完整工作流处理
"""

import os
import sys
import json
import time
import glob
import shutil
import logging
import tempfile
import subprocess
import yaml
import numpy as np
import traceback  # 添加traceback模块导入
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple, Union, Set
import re

# 导入其他必要模块
try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from molyte_cursor.src.io.inp_reader import read_inp_file, parse_inp_content
from molyte_cursor.src.utils.simulation_helpers import SimulationHelper
from molyte_cursor.src.utils.charge_modifier import replace_charges_in_lmp
from molyte_cursor.src.utils.molecule_library import MoleculeLibrary  # 导入分子库类

# 设置全局日志记录器
logger = logging.getLogger("MolyteWorkflow")

# 使用单例模式确保日志处理器只被添加一次
def setup_logger():
    """设置日志记录器，确保处理器只被添加一次"""
    # 检查是否已经有处理器
    if logger.handlers:
        return logger
        
    # 添加新的处理器
    log_handler = logging.StreamHandler()
    log_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    log_handler.setFormatter(log_formatter)
    logger.addHandler(log_handler)
    logger.setLevel(logging.INFO)
    
    # 防止日志向上传播
    logger.propagate = False
    
    # 添加任务ID到日志格式
    def add_task_id(record):
        if hasattr(record, 'task_id'):
            record.msg = f"[Task {record.task_id}] {record.msg}"
        return True
    
    logger.addFilter(add_task_id)
    
    return logger

# 初始化日志记录器
logger = setup_logger()

class MolyteWorkflow:
    """电解液计算工作流类 - 整合版"""
    
    def __init__(self, inp_file_path=None, inp_content=None, work_dir=None, output_dir=None, 
                 user_id=None, project_id=None, formulation_id=None, task_id=None, 
                 formulation_data=None, config=None):
        """初始化工作流实例
        
        Args:
            inp_file_path: 输入文件路径
            inp_content: 输入文件内容
            work_dir: 工作目录
            output_dir: 输出目录
            user_id: 用户ID
            project_id: 项目ID
            formulation_id: 配方ID
            task_id: 任务ID
            formulation_data: 配方数据，直接传入的结构化配方数据
            config: 额外配置，如auto_submit等
        """
        # 初始化日志记录器
        self.logger = setup_logger()
        
        # 初始化基本属性
        self.inp_file_path = inp_file_path
        self.inp_content = inp_content
        self.work_dir = work_dir
        self.output_dir = output_dir
        self.user_id = user_id
        self.project_id = project_id
        self.formulation_id = formulation_id
        self.task_id = task_id
        
        # 初始化结果和文件相关变量
        self.results = {}
        self.generated_files = {}
        self.pdb_file = None
        self.data_file = None
        self.in_file = None
        self.lammps_files_created = False
        
        # 初始化组件相关变量
        self.components = {
            "cations": [],
            "anions": [],
            "solvents": []
        }
        
        # 检查并记录formulation_data
        if formulation_data is not None:
            if isinstance(formulation_data, dict):
                self.logger.info(f"formulation_data 类型正确: dict")
                self.formulation_data = formulation_data
            else:
                self.logger.warning(f"formulation_data 类型错误: {type(formulation_data).__name__}, 将使用空字典")
                self.formulation_data = {
                    "name": "默认电解液配方",
                    "description": "自动生成的默认配方",
                    "cations": [],
                    "anions": [],
                    "solvents": [],
                    "parameters": {"box_size": 50.0}
                }
        else:
            self.logger.info("formulation_data 为 None，将使用空字典")
            self.formulation_data = {
                "name": "默认电解液配方",
                "description": "自动生成的默认配方",
                "cations": [],
                "anions": [],
                "solvents": [],
                "parameters": {"box_size": 50.0}
            }
        
        # 检查并记录config
        if config is not None:
            if isinstance(config, dict):
                self.logger.info(f"config 类型正确: dict")
                self.config = config
            else:
                self.logger.warning(f"config 类型错误: {type(config).__name__}, 将使用默认配置")
                self.config = {"auto_submit": True, "cpus": 64, "max_hours": 72}
        else:
            self.logger.info("config 为 None，将使用默认配置")
            self.config = {"auto_submit": True, "cpus": 64, "max_hours": 72}
        
        # 初始化工作目录
        if self.work_dir:
            os.makedirs(self.work_dir, exist_ok=True)
        
        # 初始化输出目录
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
        
        # 加载配置文件
        self.config_file = os.path.join(os.path.dirname(__file__), "../../config/paths/default_paths.yaml")
        self.logger.info(f"加载配置文件: {self.config_file}")
        try:
            with open(self.config_file, "r") as f:
                self.paths_config = yaml.safe_load(f)
            self.logger.info("成功加载配置文件")
        except Exception as e:
            self.logger.error(f"加载配置文件失败: {str(e)}")
            self.paths_config = {"paths": {}}
        
        # 初始化外部程序路径
        self._set_external_paths()
        
        # 初始化分子库
        self.molecule_library_path = self.paths_config.get("paths", {}).get("molecule_library", "molecules_library")
        self.logger.info(f"分子库初始化成功: {self.molecule_library_path}")
        
        # 检查分子库是否为空
        if os.path.exists(self.molecule_library_path) and os.listdir(self.molecule_library_path):
            self.logger.info("分子库不为空，跳过初始分子导入")
        else:
            self.logger.warning("分子库为空或不存在，后续将尝试初始化")
        
        # 记录外部程序路径
        self.logger.info("外部程序路径设置:")
        for key, value in self.paths_config.get("paths", {}).items():
            self.logger.info(f"  {key}: {value}")
        
        # 记录配置信息
        if self.config:
            self.logger.info("初始化时配置: 使用提供的config")
            self.logger.info(f"配置中的参数: {', '.join(self.config.keys())}")
        else:
            self.logger.info("初始化时配置: 未提供config")
        
        # 设置外部程序路径字典
        self.external_paths = {}
        
        # 设置是否在计算后运行模拟
        self.run_simulation = False  # 默认不运行模拟
        
        # 设置任务工作目录
        # 确保有工作空间根目录
        workspace_root = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "workspace")
        os.makedirs(workspace_root, exist_ok=True)
        
        # 设置用户工作空间
        user_workspace = os.path.join(workspace_root, "users", str(user_id) if user_id else "default")
        os.makedirs(user_workspace, exist_ok=True)
        
        # 设置任务ID目录
        task_id_dir = self.task_id
        
        # 设置工作目录
        if work_dir:
            self.work_dir = work_dir
        else:
            # 默认工作目录在workspace中
            self.work_dir = os.path.join(user_workspace, task_id_dir)
        
        # 设置输出目录
        self.output_dir = output_dir if output_dir else self.work_dir
        
        # 创建工作目录和输出目录
        os.makedirs(self.work_dir, exist_ok=True)
        if self.output_dir != self.work_dir:
            os.makedirs(self.output_dir, exist_ok=True)
            
        self.logger.info(f"工作目录: {self.work_dir}")
        self.logger.info(f"输出目录: {self.output_dir}")
        
        # 存储输入文件路径和内容
        self.inp_file_path = inp_file_path
        self.inp_content = inp_content
        
        # 初始化结果和文件路径
        # 如果提供了配置，直接使用它
        self.config = config if config is not None else {}
        self.logger.info(f"初始化时配置: {'使用提供的config' if config else '使用空配置'}")
        if config:
            self.logger.info(f"配置中的参数: {', '.join(config.keys())}")
            if 'box_size' in config:
                self.logger.info(f"配置中的box_size: {config.get('box_size')} ({type(config.get('box_size')).__name__})")
            if 'auto_submit' in config:
                self.logger.info(f"配置中的auto_submit: {config.get('auto_submit')}")
        
        self.molecule_files = {}  # 分子文件
        self.generated_files = {}  # 生成的文件
        self.results = {}  # 计算结果
    
    def _set_external_paths(self):
        """设置外部程序和路径"""
        paths = self.paths_config.get('paths', {})
        
        # 设置外部程序路径
        self.ligpargen_path = paths.get('ligpargen')
        self.packmol_path = paths.get('packmol')
        self.ltemplify_path = paths.get('ltemplify')
        self.moltemplate_path = paths.get('moltemplate')
        self.resp_path = paths.get('resp')
        
        # 设置数据目录
        self.common_data_path = paths.get('common_data')
        self.charge_save_path = paths.get('charge_save')
        self.molecule_workspace = paths.get('molecule_workspace')
        
        # 初始分子目录
        molyte_base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.initial_salts_path = paths.get('initial_salts', os.path.join(molyte_base_path, 'data', 'initial_salts'))
        self.initial_solvent_path = paths.get('initial_solvent', os.path.join(molyte_base_path, 'data', 'initial_solvent'))
        
        # 设置分子库路径
        self.molecules_library_path = paths.get('molecules_library', os.path.join(molyte_base_path, 'molecules_library'))
        
        # 确保分子库目录结构存在
        try:
            # 创建分子库主目录
            os.makedirs(self.molecules_library_path, exist_ok=True)
            
            # 创建必要的子目录结构
            metadata_dir = os.path.join(self.molecules_library_path, "metadata")
            ions_dir = os.path.join(self.molecules_library_path, "ions")
            cations_dir = os.path.join(ions_dir, "cations")
            anions_dir = os.path.join(ions_dir, "anions")
            solvents_dir = os.path.join(self.molecules_library_path, "solvents")
            
            os.makedirs(metadata_dir, exist_ok=True)
            os.makedirs(cations_dir, exist_ok=True)
            os.makedirs(anions_dir, exist_ok=True)
            os.makedirs(solvents_dir, exist_ok=True)
            
            # 确保元数据文件存在
            molecules_json = os.path.join(metadata_dir, "molecules.json")
            smiles_index_json = os.path.join(metadata_dir, "smiles_index.json")
            
            if not os.path.exists(molecules_json):
                with open(molecules_json, 'w') as f:
                    json.dump({"ions": {"cations": [], "anions": []}, "solvents": []}, f, indent=2)
                self.logger.info(f"创建了新的分子库元数据文件: {molecules_json}")
                
            if not os.path.exists(smiles_index_json):
                with open(smiles_index_json, 'w') as f:
                    json.dump({}, f, indent=2)
                self.logger.info(f"创建了新的SMILES索引文件: {smiles_index_json}")
            
            # 导入MoleculeLibrary模块
            from src.utils.molecule_library import MoleculeLibrary
            
            # 初始化分子库
            self.molecule_library = MoleculeLibrary(self.molecules_library_path)
            self.logger.info(f"分子库初始化成功: {self.molecules_library_path}")
            self.use_molecule_library = True
            
            # 导入初始分子到分子库(如果需要)
            self._import_initial_molecules_to_library()
        except ImportError as e:
            self.logger.warning(f"导入MoleculeLibrary模块失败: {str(e)}，将使用传统文件处理")
            self.molecule_library = None
            self.use_molecule_library = False
        except Exception as e:
            self.logger.warning(f"分子库初始化失败: {str(e)}，将使用传统文件处理")
            self.molecule_library = None
            self.use_molecule_library = False
        
        # 是否使用RESP电荷
        self.use_resp_charges = paths.get('use_resp_charges', True)
        
        # 加载分子库文件
        self.exist_molecule_path = os.path.join(self.common_data_path, 'exist_molecule.json') if self.common_data_path else None
        
        # 记录路径信息
        self.logger.info("外部程序路径设置:")
        self.logger.info(f"  ligpargen: {self.ligpargen_path}")
        self.logger.info(f"  packmol: {self.packmol_path}")
        self.logger.info(f"  ltemplify: {self.ltemplify_path}")
        self.logger.info(f"  moltemplate: {self.moltemplate_path}")
        self.logger.info(f"  resp: {self.resp_path}")
        self.logger.info(f"  通用数据路径: {self.common_data_path}")
        self.logger.info(f"  电荷保存路径: {self.charge_save_path}")
        self.logger.info(f"  分子库文件: {self.exist_molecule_path}")
        self.logger.info(f"  初始盐路径: {self.initial_salts_path}")
        self.logger.info(f"  初始溶剂路径: {self.initial_solvent_path}")
        self.logger.info(f"  分子库路径: {self.molecules_library_path}")
        self.logger.info(f"  使用分子库: {self.use_molecule_library}")
    
    def _import_initial_molecules_to_library(self):
        """导入初始分子到分子库
        
        从初始盐和溶剂目录导入分子到分子库，仅在分子库为空时执行
        """
        if not hasattr(self, 'molecule_library') or not self.molecule_library:
            return
            
        # 检查分子库是否为空
        molecules = self.molecule_library.list_molecules()
        if molecules.get("ions", {}).get("cations") or molecules.get("ions", {}).get("anions") or molecules.get("solvents"):
            self.logger.info("分子库不为空，跳过初始分子导入")
            return
            
        self.logger.info("检测到空分子库，开始导入初始分子")
        
        # 导入初始盐
        if os.path.exists(self.initial_salts_path):
            cation_names = []
            anion_names = []
            
            # 为了简单起见，假设所有文件名都是分子名称
            for filename in os.listdir(self.initial_salts_path):
                if filename.endswith('.pdb'):
                    base_name = os.path.splitext(filename)[0]
                    # 简单规则：带+的是阳离子，带-的是阴离子
                    if '+' in base_name:
                        cation_names.append(base_name)
                    elif '-' in base_name:
                        anion_names.append(base_name)
            
            # 导入阳离子
            for name in cation_names:
                files = {}
                for ext in ['pdb', 'lt', 'lmp']:
                    file_path = os.path.join(self.initial_salts_path, f"{name}.{ext}")
                    if os.path.exists(file_path):
                        files[ext] = file_path
                
                if files:
                    try:
                        self.molecule_library.add_molecule(name, files, "cations")
                        self.logger.info(f"已导入初始阳离子: {name}")
                    except Exception as e:
                        self.logger.warning(f"导入初始阳离子 {name} 失败: {str(e)}")
            
            # 导入阴离子
            for name in anion_names:
                files = {}
                for ext in ['pdb', 'lt', 'lmp']:
                    file_path = os.path.join(self.initial_salts_path, f"{name}.{ext}")
                    if os.path.exists(file_path):
                        files[ext] = file_path
                
                if files:
                    try:
                        self.molecule_library.add_molecule(name, files, "anions")
                        self.logger.info(f"已导入初始阴离子: {name}")
                    except Exception as e:
                        self.logger.warning(f"导入初始阴离子 {name} 失败: {str(e)}")
        
        # 导入初始溶剂
        if os.path.exists(self.initial_solvent_path):
            # 获取初始溶剂目录下的所有.pdb文件
            solvent_files = [f for f in os.listdir(self.initial_solvent_path) if f.endswith('.pdb')]
            
            for pdb_file in solvent_files:
                name = os.path.splitext(pdb_file)[0]
                files = {}
                
                # 收集所有相关文件
                for ext in ['pdb', 'lt', 'lmp', 'chg']:
                    file_path = os.path.join(self.initial_solvent_path, f"{name}.{ext}")
                    if os.path.exists(file_path):
                        files[ext] = file_path
                
                # 查找或创建SMILE文件
                smile_file = os.path.join(self.initial_solvent_path, f"{name}.smile")
                smile = None
                
                if os.path.exists(smile_file):
                    with open(smile_file, 'r') as f:
                        smile = f.read().strip()
                
                if files and smile:
                    try:
                        self.molecule_library.add_molecule(name, files, "solvents", smile)
                        self.logger.info(f"已导入初始溶剂: {name}, SMILE: {smile}")
                    except Exception as e:
                        self.logger.warning(f"导入初始溶剂 {name} 失败: {str(e)}")
                elif files:
                    self.logger.warning(f"溶剂 {name} 缺少SMILE，无法导入到分子库")
        
        self.logger.info("初始分子导入完成")
    
    @classmethod
    def from_formulation(cls, formulation_data=None, work_dir=None, output_dir=None,
                        user_id=None, project_id=None, formulation_id=None,
                        task_id=None, config=None):
        """从配方数据创建工作流实例
        
        Args:
            formulation_data: 配方数据，包含cations, anions, solvents等
            work_dir: 工作目录
            output_dir: 输出目录
            user_id: 用户ID
            project_id: 项目ID
            formulation_id: 配方ID
            task_id: 任务ID
            config: 额外的配置参数，如auto_submit等
        
        Returns:
            MolyteWorkflow: 工作流实例
        """
        # 创建日志记录器
        logger = setup_logger()
        
        # 记录传入的参数
        logger.info(f"开始创建工作流: formulation_id={formulation_id}, task_id={task_id}")
        
        # 参数检查 - 记录formulation_data
        logger.info(f"【参数调试】接收到的完整formulation_data数据: {formulation_data}")
        
        # 确保formulation_data不为None
        if formulation_data is None:
            logger.warning("formulation_data为None，初始化为默认配置")
            formulation_data = {
                "name": "默认电解液配方",
                "description": "自动生成的默认配方",
                "parameters": {"box_size": 50.0},
                "cations": [],
                "anions": [],
                "solvents": []
            }
        
        # 从formulation_data中提取关键参数
        box_size = None
        if isinstance(formulation_data, dict) and 'parameters' in formulation_data and 'box_size' in formulation_data['parameters']:
            box_size = formulation_data['parameters']['box_size']
            try:
                box_size = float(box_size)
                logger.info(f"【关键参数】从formulation中提取的box_size: {box_size}")
            except (ValueError, TypeError) as e:
                logger.warning(f"box_size转换为float失败: {e}，使用默认值")
                box_size = 50.0
        
        # 记录配方的基本信息
        if isinstance(formulation_data, dict):
            logger.info(f"成功从数据库构建配方数据: {formulation_data.get('name', '未命名配方')}")
            
            # 记录组件数量
            logger.info(f"阳离子数量: {len(formulation_data.get('cations', []))}")
            logger.info(f"阴离子数量: {len(formulation_data.get('anions', []))}")
            logger.info(f"溶剂数量: {len(formulation_data.get('solvents', []))}")
            
            # 提取模拟参数
            simulation_params = formulation_data.get('parameters', {})
            logger.info(f"【参数调试】准备添加模拟参数: {simulation_params}")
            
            # 确保box_size是float类型
            if 'box_size' in simulation_params:
                try:
                    simulation_params['box_size'] = float(simulation_params['box_size'])
                    logger.info(f"【关键参数】已将box_size转换为float: {simulation_params['box_size']}")
                except (ValueError, TypeError) as e:
                    logger.warning(f"转换box_size失败: {e}，使用默认值")
                    simulation_params['box_size'] = 50.0
        
        # 应用config参数（如果提供）
        if config is None:
            # 创建默认配置
            config = {
                "auto_submit": True,
                "cpus": 64,
                "max_hours": 72
            }
        elif not isinstance(config, dict):
            logger.warning(f"config参数类型错误: {type(config).__name__}，将使用默认配置")
            config = {
                "auto_submit": True,
                "cpus": 64,
                "max_hours": 72
            }
        else:
            # 确保config有必要的键
            if 'auto_submit' not in config:
                config['auto_submit'] = True
        
        logger.info(f"最终配置参数: {config}")
        
        # 创建工作流实例
        workflow = cls(
            inp_file_path=None,
            inp_content=None,
            work_dir=work_dir,
            output_dir=output_dir,
            user_id=user_id,
            project_id=project_id,
            formulation_id=formulation_id,
            task_id=task_id,
            formulation_data=formulation_data,
            config=config
        )
        
        return workflow

    def _execute_packmol(self, packmol_input_file, output_pdb, output_dir):
        """执行Packmol命令，生成系统PDB文件
        
        Args:
            packmol_input_file: Packmol输入文件路径
            output_pdb: 输出PDB文件路径
            output_dir: 输出目录
            
        Returns:
            Dict: 包含执行结果的字典
        """
        self.logger.info(f"开始执行Packmol: {packmol_input_file}")
        
        # 执行Packmol
        packmol_path = self.external_paths.get('packmol', 'packmol')
        
        # 为了避免输入重定向问题，使用临时文件作为输入
        packmol_cmd = f"{packmol_path} < {packmol_input_file}"
        self.logger.info(f"执行Packmol命令: {packmol_cmd}")
        
        # 使用subprocess执行外部命令，捕获输出和错误
        import subprocess
        process = subprocess.run(
            packmol_cmd, 
            shell=True, 
            cwd=output_dir,
            capture_output=True,
            text=True
        )
        
        # 记录标准输出和错误输出
        if process.stdout:
            self.logger.info(f"Packmol输出:\n{process.stdout}")
        if process.stderr:
            self.logger.warning(f"Packmol错误:\n{process.stderr}")
        
        # 检查Packmol执行结果
        if process.returncode != 0:
            self.logger.error(f"Packmol执行失败，返回码: {process.returncode}")
            return {'error': f"Packmol执行失败: {process.stderr}"}
        
        # 检查输出文件是否存在且大小合理
        if not os.path.exists(output_pdb):
            self.logger.error(f"Packmol未能生成输出文件: {output_pdb}")
            return {'error': "Packmol未能生成输出文件"}
        
        file_size = os.path.getsize(output_pdb)
        if file_size < 100:  # 文件太小，可能是空的或只有头部
            self.logger.error(f"Packmol生成的文件太小 ({file_size} bytes): {output_pdb}")
            return {'error': "Packmol生成的文件不完整"}
        
        self.logger.info(f"Packmol成功构建系统，输出文件: {output_pdb}, 大小: {file_size} bytes")
        return {'output_pdb': output_pdb}
    
    def run(self):
        """运行工作流
        
        Returns:
            Dict: 运行结果
        """
        self.logger.info("开始运行Molyte工作流")
        
        # 初始化结果字典
        if not hasattr(self, 'results') or self.results is None:
            self.results = {}
        
        self.results.update({
            "status": "running",
            "output_files": {}
        })
        
        try:
            # 解析输入文件
            self.logger.info("开始解析输入文件")
            try:
                parse_result = self.parse_input()
                self.logger.info(f"输入文件解析完成")
            except Exception as parse_error:
                self.logger.error(f"解析输入文件失败: {str(parse_error)}")
                self.logger.error(f"详细错误信息: {traceback.format_exc()}")
                raise RuntimeError(f"解析输入文件失败: {str(parse_error)}")
            
            # 创建输入文件（如果需要）
            if not hasattr(self, 'lammps_files_created') or not self.lammps_files_created:
                self.logger.info("开始创建LAMMPS输入文件")
                self.create_lammps_files()
                self.logger.info("LAMMPS输入文件创建完成")
            else:
                self.logger.info("LAMMPS输入文件已存在，跳过创建")
                
            # 如果配置中设置了auto_submit=True，则自动提交作业
            auto_submit = False
            if hasattr(self, 'config') and isinstance(self.config, dict):
                auto_submit = self.config.get("auto_submit", False)
                
            if auto_submit:
                self.logger.info("自动提交作业设置已启用")
                slurm_job_id = self.submit_job()
                if slurm_job_id:
                    self.results['slurm_job_id'] = slurm_job_id
                    self.logger.info(f"作业已提交，SLURM作业ID: {slurm_job_id}")
                else:
                    self.logger.warning("提交作业未返回有效的作业ID")
            else:
                self.logger.info("自动提交作业设置已禁用，需要手动提交")
            
            # 更新状态
            self.results['status'] = 'success'
            
            # 添加系统PDB文件路径
            if hasattr(self, 'pdb_file') and self.pdb_file and os.path.exists(self.pdb_file):
                self.results['system_pdb'] = self.pdb_file
                
            # 添加其他输出文件路径
            if hasattr(self, 'generated_files') and self.generated_files:
                self.results['output_files'] = self.generated_files
                
            self.logger.info("工作流成功完成")
            return self.results
            
        except Exception as e:
            self.logger.error(f"工作流执行失败: {str(e)}")
            import traceback
            self.logger.error(f"详细错误信息: {traceback.format_exc()}")
            
            # 更新状态
            if not hasattr(self, 'results') or self.results is None:
                self.results = {}
                
            self.results.update({
                "status": "failed",
                "error": str(e)
            })
            
            return self.results
    
    def collect_results(self) -> Dict[str, Any]:
        """收集计算结果
        
        收集并返回电解质计算的结果数据
        
        Returns:
            Dict[str, Any]: 计算结果字典
        """
        self.logger.info("收集计算结果")
        
        # 如果结果为空，返回默认值
        if not hasattr(self, 'results') or not self.results:
            self.logger.warning("没有找到计算结果，返回默认值")
            return {
                "status": "completed",
                "conductivity": 0.0,
                "diffusion_coefficients": {
                    "cation": 0.0,
                    "anion": 0.0
                },
                "density": 0.0,
                "viscosity": 0.0,
                "timestamp": datetime.now().isoformat()
            }
        
        # 返回已存储的结果
        return self.results
    
    def parse_input(self):
        """解析输入文件"""
        self.logger.info("开始解析输入文件")
        self.logger.info("调试信息 - 对象属性:")
        self.logger.info(f"has formulation_data: {hasattr(self, 'formulation_data')}")
        
        # 确保结果字典已初始化
        if not hasattr(self, 'results') or self.results is None:
            self.results = {}
        
        # 确保generated_files已初始化
        if not hasattr(self, 'generated_files') or self.generated_files is None:
            self.generated_files = {}
        
        # 增加formulation_data类型调试
        if hasattr(self, 'formulation_data'):
            self.logger.info(f"formulation_data 值类型: {type(self.formulation_data)}")
            
            # 防御性检查：确保self.formulation_data不为None
            if self.formulation_data is None:
                self.logger.warning("发现formulation_data为None，将初始化为空字典")
                self.formulation_data = {
                    "name": "默认电解液配方",
                    "description": "自动生成的默认配方",
                    "cations": [],
                    "anions": [],
                    "solvents": [],
                    "parameters": {"box_size": 50.0}
                }
        else:
            self.logger.warning("对象没有formulation_data属性，将创建默认值")
            self.formulation_data = {
                "name": "默认电解液配方",
                "description": "自动生成的默认配方",
                "cations": [],
                "anions": [],
                "solvents": [],
                "parameters": {"box_size": 50.0}
            }
        
        # 确保components已初始化
        if not hasattr(self, 'components') or self.components is None:
            self.components = {
                "cations": [],
                "anions": [],
                "solvents": []
            }
        
        # 从formulation_data中提取组件信息（只有在确保它是有效的字典后）
        if isinstance(self.formulation_data, dict):
            if 'cations' in self.formulation_data:
                self.logger.info(f"formulation_data 中 cations 长度: {len(self.formulation_data['cations'])}")
            else:
                self.logger.info("formulation_data 中没有 cations 键")
                # 确保存在cations字段，即使为空
                self.formulation_data['cations'] = []
            
            if 'anions' in self.formulation_data:
                self.logger.info(f"formulation_data 中 anions 长度: {len(self.formulation_data['anions'])}")
            else:
                self.logger.info("formulation_data 中没有 anions 键")
                # 确保存在anions字段，即使为空
                self.formulation_data['anions'] = []
            
            if 'solvents' in self.formulation_data:
                self.logger.info(f"formulation_data 中 solvents 长度: {len(self.formulation_data['solvents'])}")
            else:
                self.logger.info("formulation_data 中没有 solvents 键")
                # 确保存在solvents字段，即使为空
                self.formulation_data['solvents'] = []
            
            # 确保parameters也存在
            if 'parameters' not in self.formulation_data:
                self.logger.info("formulation_data 中没有 parameters 键，添加默认参数")
                self.formulation_data['parameters'] = {"box_size": 50.0}
        else:
            self.logger.warning(f"formulation_data 不是字典类型，重新初始化: {type(self.formulation_data)}")
            self.formulation_data = {
                "name": "默认电解液配方",
                "description": "自动生成的默认配方",
                "cations": [],
                "anions": [],
                "solvents": [],
                "parameters": {"box_size": 50.0}
            }
        
        # 安全地更新结果字典
        self.results['output_files'] = self.generated_files
        
        return self.results
    
    def _parse_key_value_format(self, file_path):
        """解析键值对格式的INP文件
        
        Args:
            file_path: INP文件路径
            
        Returns:
            dict: 解析后的数据字典
        """
        self.logger.info(f"解析键值对格式文件: {file_path}")
        data = {}
        inside_block = False
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 跳过空行和注释行
                if not line or line.startswith('#'):
                    continue
                    
                # 处理START和END标记
                if line.upper() == 'START':
                    inside_block = True
                    continue
                elif line.upper() == 'END':
                    inside_block = False
                    continue
                
                # 只在START和END之间处理内容
                if inside_block:
                    # 分割键值对
                    parts = line.split(None, 1)
                    if len(parts) == 2:
                        key, value = parts
                        data[key] = value
        
        # 提取分子信息
        result = self._extract_molecules_from_key_value(data)
        
        return result
        
    def _extract_molecules_from_key_value(self, data):
        """从键值对数据中提取分子信息
        
        Args:
            data: 键值对数据字典
            
        Returns:
            dict: 包含结构化分子信息的数据字典
        """
        # 初始化配置字典
        result = {
            'name': data.get('name', '未命名配方'),
            'temperature': float(data.get('T', 298.15)),
            'box_size': float(data.get('box_size', 50.0)),  # 使用小写的box_size字段名
            'concentration': float(data.get('concentration', 1.0)),
            'time_step': float(data.get('time_step', 1.0)),
            'equilibration_steps': int(data.get('equilibration_steps', 1000000)),
            'production_steps': int(data.get('production_steps', 2000000)),
            'cutoff': float(data.get('cutoff', 12.0)),
            'pressure': float(data.get('pressure', 1.0)),
            'cations': [],
            'anions': [],
            'solvents': []
        }
        
        # 获取盒子大小并确保有效值
        # 支持box_size和Box_size两种情况，优先使用小写版本
        orig_box_size_lower = data.get('box_size')
        orig_box_size_upper = data.get('Box_size')
        self.logger.info(f"【跟踪box_size】原始数据中的box_size(小写): {orig_box_size_lower}")
        self.logger.info(f"【跟踪box_size】原始数据中的Box_size(大写): {orig_box_size_upper}")
        
        box_size = max(float(data.get('box_size', data.get('Box_size', 50.0))), 10.0)  # 确保至少是10Å
        
        # 更新结果字典中的box_size
        result['box_size'] = box_size
        
        # 记录处理的盒子大小，便于调试
        self.logger.info(f"【跟踪box_size】_extract_molecules_from_key_value方法中最终使用的box_size: {box_size}")
        self.logger.info(f"盒子大小: {box_size} Å")
        
        # 计算盒子体积 (Å³)
        box_volume_A3 = box_size ** 3
        
        # 防御性检查：确保体积不为零
        if box_volume_A3 <= 0:
            self.logger.warning(f"盒子体积计算为零或负值: {box_volume_A3}，使用默认值代替")
            box_volume_A3 = 50.0 ** 3  # 使用默认50Å盒子的体积
        
        # 转换为升 (L), 1 Å³ = 10^-27 L
        box_volume_L = box_volume_A3 * 1e-27
        
        # 记录体积信息，便于调试
        self.logger.info(f"盒子体积: {box_volume_L:.6f} L")
        
        # 获取浓度 (mol/L)
        concentration = float(data.get('concentration', 1.0))
        
        # 阿伏伽德罗常数 (个/mol)
        AVOGADRO = 6.02214076e23
        
        # 计算基本分子数量（阳离子的数量，根据浓度和体积确定）
        # 分子数 = 浓度(mol/L) × 体积(L) × 阿伏伽德罗常数(个/mol)
        # 使用try-except防止数值溢出
        try:
            base_molecules = int(concentration * box_volume_L * AVOGADRO)
        except (OverflowError, ValueError):
            self.logger.warning("计算分子数量时发生溢出，使用默认值代替")
            base_molecules = 16  # 默认值
        
        # 如果计算出的分子数量不合理（过大或过小），设置合理的默认值
        if base_molecules <= 0:
            self.logger.warning(f"计算的分子数量不合理: {base_molecules}，使用默认值16")
            base_molecules = 16  # 默认每个离子有16个分子
        elif base_molecules > 100:
            # 限制最大分子数以避免计算资源问题
            old_value = base_molecules
            base_molecules = 100
            self.logger.warning(f"计算的分子数量过大: {old_value}，限制为最大值: {base_molecules}")
        
        # 记录计算的分子数量，便于调试
        self.logger.info(f"基础分子数量计算结果: {base_molecules}")
        
        # 从键值对中提取阳离子信息
        i = 1
        cation_total = 0  # 跟踪阳离子总数
        
        # 支持web格式的cation字段
        if any(k.startswith('cation') for k in data.keys()):
            for k, v in data.items():
                if k.startswith('cation'):
                    parts = v.split(',')
                    if len(parts) >= 2:
                        cation_name = parts[0].strip()
                        try:
                            # 优先使用前端提交的数量
                            count = int(parts[1].strip())
                        except ValueError:
                            # 如果解析失败，使用基于计算的默认值
                            count = base_molecules
                            self.logger.info(f"阳离子 {cation_name} 未指定有效数量，使用计算值: {count}")
                        
                        # 安全检查，避免数值过大
                        count = min(count, 100)  # 最多100个阳离子分子
                        
                        result['cations'].append({
                            'name': cation_name,
                            'ratio': 1.0,
                            'count': count,
                            'type': 'cation'
                        })
                        cation_total += count
                        self.logger.info(f"添加阳离子 {cation_name}，数量: {count}")
        else:
            # 原有的Cation1_name格式支持
            while True:
                cation_name_key = f"Cation{i}_name"
                cation_count_key = f"Cation{i}_count"  # 尝试先查找count键
                cation_ratio_key = f"Cation{i}_ratio"
                
                if cation_name_key not in data:
                    break
                    
                cation_name = data.get(cation_name_key)
                
                # 优先使用直接提供的count
                if cation_count_key in data:
                    try:
                        count = int(data.get(cation_count_key))
                    except ValueError:
                        # 如果解析失败，使用基于比例的计算方式
                        cation_ratio = float(data.get(cation_ratio_key, 1.0))
                        count = max(1, int(base_molecules * cation_ratio))
                        self.logger.info(f"阳离子 {cation_name} count值无效，使用比例计算: {count}")
                else:
                    # 使用比例计算
                    cation_ratio = float(data.get(cation_ratio_key, 1.0))
                    count = max(1, int(base_molecules * cation_ratio))
                    self.logger.info(f"阳离子 {cation_name} 使用比例计算数量: {count}")
                
                # 设置最大值限制，避免数值过大
                count = min(count, 100)  # 设置上限为100个分子
                
                result['cations'].append({
                    'name': cation_name,
                    'ratio': float(data.get(cation_ratio_key, 1.0)),
                    'count': count,
                    'type': 'cation'
                })
                
                cation_total += count
                self.logger.info(f"添加阳离子 {cation_name}，数量: {count}")
                i += 1
        
        # 记录阳离子信息
        self.logger.info(f"阳离子总数: {cation_total}")
            
        # 从键值对中提取阴离子信息
        i = 1
        anion_total = 0  # 跟踪阴离子总数
        
        # 支持web格式的anion字段
        if any(k.startswith('anion') for k in data.keys()):
            for k, v in data.items():
                if k.startswith('anion'):
                    parts = v.split(',')
                    if len(parts) >= 2:
                        anion_name = parts[0].strip()
                        try:
                            # 优先使用前端提交的数量
                            count = int(parts[1].strip())
                        except ValueError:
                            # 如果解析失败，使用基于计算的默认值
                            count = base_molecules
                            self.logger.info(f"阴离子 {anion_name} 未指定有效数量，使用计算值: {count}")
                        
                        # 安全检查，避免数值过大
                        count = min(count, 100)  # 最多100个阴离子分子
                        
                        result['anions'].append({
                            'name': anion_name,
                            'ratio': 1.0,
                            'count': count,
                            'type': 'anion'
                        })
                        anion_total += count
                        self.logger.info(f"添加阴离子 {anion_name}，数量: {count}")
        else:
            # 原有的Anion1_name格式支持
            while True:
                anion_name_key = f"Anion{i}_name"
                anion_count_key = f"Anion{i}_count"  # 尝试先查找count键
                anion_ratio_key = f"Anion{i}_ratio"
                
                if anion_name_key not in data:
                    break
                    
                anion_name = data.get(anion_name_key)
                
                # 优先使用直接提供的count
                if anion_count_key in data:
                    try:
                        count = int(data.get(anion_count_key))
                    except ValueError:
                        # 如果解析失败，使用基于比例的计算方式
                        anion_ratio = float(data.get(anion_ratio_key, 1.0))
                        count = max(1, int(base_molecules * anion_ratio))
                        self.logger.info(f"阴离子 {anion_name} count值无效，使用比例计算: {count}")
                else:
                    # 使用比例计算
                    anion_ratio = float(data.get(anion_ratio_key, 1.0))
                    count = max(1, int(base_molecules * anion_ratio))
                    self.logger.info(f"阴离子 {anion_name} 使用比例计算数量: {count}")
                
                # 设置最大值限制，避免数值过大
                count = min(count, 100)  # 设置上限为100个分子
                
                result['anions'].append({
                    'name': anion_name,
                    'ratio': float(data.get(anion_ratio_key, 1.0)),
                    'count': count,
                    'type': 'anion'
                })
                
                anion_total += count
                self.logger.info(f"添加阴离子 {anion_name}，数量: {count}")
                i += 1
        
        # 记录阴离子信息
        self.logger.info(f"阴离子总数: {anion_total}")
            
        # 从键值对中提取溶剂信息
        i = 1
        # 溶剂数量通常是离子的3-4倍
        solvent_base = base_molecules * 3
        solvent_total = 0  # 跟踪溶剂总数
        
        # 支持web格式的solvent字段
        if any(k.startswith('solvent') for k in data.keys()):
            for k, v in data.items():
                if k.startswith('solvent'):
                    parts = v.split(',')
                    if len(parts) >= 2:
                        solvent_name = parts[0].strip()
                        # 溶剂可能还包含SMILE字符串
                        solvent_smile = ""
                        
                        # 处理各种格式：name,smile,count 或 name,count
                        if len(parts) > 2:
                            # 假设格式为 name,smile,count
                            solvent_smile = parts[1].strip()
                            try:
                                count = int(parts[2].strip())
                            except (ValueError, IndexError):
                                # 如果解析失败，使用默认值
                                count = solvent_base
                                self.logger.info(f"溶剂 {solvent_name} 未指定有效数量，使用计算值: {count}")
                        else:
                            # 假设格式为 name,count
                            try:
                                count = int(parts[1].strip())
                            except ValueError:
                                # 如果解析失败，使用默认值
                                count = solvent_base
                                self.logger.info(f"溶剂 {solvent_name} 未指定有效数量，使用计算值: {count}")
                        
                        # 安全检查，避免数值过大
                        count = min(count, 300)  # 最多300个溶剂分子
                        
                        result['solvents'].append({
                            'name': solvent_name,
                            'smile': solvent_smile,
                            'ratio': 1.0,
                            'count': count,
                            'type': 'solvent'
                        })
                        solvent_total += count
                        self.logger.info(f"添加溶剂 {solvent_name}，数量: {count}")
        else:
            # 原有的Sol1_name格式支持
            while True:
                solvent_name_key = f"Sol{i}_name"
                solvent_smile_key = f"Sol{i}_smile"
                solvent_count_key = f"Sol{i}_count"  # 尝试先查找count键
                solvent_ratio_key = f"Sol{i}_ratio"
                
                if solvent_name_key not in data:
                    break
                    
                solvent_name = data.get(solvent_name_key)
                solvent_smile = data.get(solvent_smile_key, "")
                
                # 优先使用直接提供的count
                if solvent_count_key in data:
                    try:
                        count = int(data.get(solvent_count_key))
                    except ValueError:
                        # 如果解析失败，使用基于比例的计算方式
                        solvent_ratio = float(data.get(solvent_ratio_key, 1.0))
                        count = max(4, int(solvent_base * solvent_ratio))
                        self.logger.info(f"溶剂 {solvent_name} count值无效，使用比例计算: {count}")
                else:
                    # 使用比例计算
                    solvent_ratio = float(data.get(solvent_ratio_key, 1.0))
                    count = max(4, int(solvent_base * solvent_ratio))
                    self.logger.info(f"溶剂 {solvent_name} 使用比例计算数量: {count}")
                
                # 设置最大值限制，避免数值过大
                count = min(count, 300)  # 设置上限为300个溶剂分子
                
                result['solvents'].append({
                    'name': solvent_name,
                    'smile': solvent_smile,
                    'ratio': float(data.get(solvent_ratio_key, 1.0)),
                    'count': count,
                    'type': 'solvent'
                })
                
                solvent_total += count
                self.logger.info(f"添加溶剂 {solvent_name}，数量: {count}, SMILE: {solvent_smile}")
                i += 1
        
        # 记录溶剂信息
        self.logger.info(f"溶剂总数: {solvent_total}")
        
        # 处理计算类型
        if 'calculation_types' in data:
            calculation_types = data['calculation_types']
            if isinstance(calculation_types, str):
                calculation_types = [ct.strip() for ct in calculation_types.split(',')]
            result['calculation_types'] = calculation_types
        else:
            result['calculation_types'] = []
            
        return result
    
    def generate_input_files(self):
        """生成计算所需的输入文件
        
        生成packmol和lammps输入文件，并处理所有分子文件
        """
        self.logger.info("开始生成输入文件")
        
        # 设置使用分子库标志
        self.use_molecule_library = hasattr(self, 'molecule_library') and self.molecule_library is not None
        self.logger.info(f"使用分子库: {self.use_molecule_library}")
        
        # 创建分子文件目录
        molecule_dir = os.path.join(self.work_dir, "molecules")
        os.makedirs(molecule_dir, exist_ok=True)
        
        # 处理所有分子
        self.logger.info("处理分子文件")
        molecules_info = []
        
        # 处理阳离子
        for cation in self.formulation_data.get("cations", []):
            molecule_name = cation.get("name")
            self.logger.info(f"处理阳离子: {molecule_name}")
            files = self._generate_molecule_files(cation, molecule_dir)
            if files:
                sanitized_name = molecule_name.replace('+', '').replace('-', '')
                cation["sanitized_name"] = sanitized_name
                cation["files"] = files
                molecules_info.append(cation)
                self.logger.info(f"阳离子 {molecule_name} 文件处理完成")
            else:
                self.logger.error(f"阳离子 {molecule_name} 文件生成失败")
                
        # 处理阴离子
        for anion in self.formulation_data.get("anions", []):
            molecule_name = anion.get("name")
            self.logger.info(f"处理阴离子: {molecule_name}")
            files = self._generate_molecule_files(anion, molecule_dir)
            if files:
                sanitized_name = molecule_name.replace('+', '').replace('-', '')
                anion["sanitized_name"] = sanitized_name
                anion["files"] = files
                molecules_info.append(anion)
                self.logger.info(f"阴离子 {molecule_name} 文件处理完成")
            else:
                self.logger.error(f"阴离子 {molecule_name} 文件生成失败")
                
        # 处理溶剂
        for solvent in self.formulation_data.get("solvents", []):
            molecule_name = solvent.get("name")
            self.logger.info(f"处理溶剂: {molecule_name}")
            files = self._generate_molecule_files(solvent, molecule_dir)
            if files:
                sanitized_name = molecule_name.replace('+', '').replace('-', '')
                solvent["sanitized_name"] = sanitized_name
                solvent["files"] = files
                molecules_info.append(solvent)
                self.logger.info(f"溶剂 {molecule_name} 文件处理完成")
            else:
                self.logger.error(f"溶剂 {molecule_name} 文件生成失败")
        
        # 如果没有任何分子处理成功，抛出异常
        if not molecules_info:
            raise RuntimeError("未能成功处理任何分子文件")
            
        # 更新分子信息
        self.molecule_files = {molecule.get("name"): molecule.get("files", {}) for molecule in molecules_info}
        self.logger.info(f"成功处理 {len(self.molecule_files)} 个分子文件")
            
        # 构建完整系统
        output_dir = os.path.join(self.work_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger.info("开始构建系统")
        system_files = self._build_system_with_packmol(molecules_info, output_dir)
        
        if not system_files or not system_files.get("pdb"):
            raise RuntimeError("系统构建失败")
            
        self.generated_files = {
            "system_pdb": system_files.get("pdb"),
            "system_lt": system_files.get("lt"),
            "lammps_input": system_files.get("in"),
            "lammps_data": system_files.get("data")
        }
        
        self.logger.info("输入文件生成完成")
        self.logger.info(f"系统PDB文件: {self.generated_files['system_pdb']}")
        self.logger.info(f"系统LT文件: {self.generated_files['system_lt']}")
        self.logger.info(f"LAMMPS输入文件: {self.generated_files['lammps_input']}")
        self.logger.info(f"LAMMPS数据文件: {self.generated_files['lammps_data']}")
        
        return self.generated_files
    
    def _generate_molecule_files(self, molecule_info, molecule_dir):
        """生成或复制分子文件"""
        molecule_name = molecule_info.get('name', '')
        molecule_smile = molecule_info.get('smile', '')
        
        # 确定分子类别
        category = None
        for cat in self.formulation_data.get("cations", []):
            if cat.get("name") == molecule_name:
                category = "cations"
                molecule_info['type'] = 'cation'
                break
        
        if not category:
            for an in self.formulation_data.get("anions", []):
                if an.get("name") == molecule_name:
                    category = "anions"
                    molecule_info['type'] = 'anion'
                    break
        
        if not category:
            category = "solvents"
            molecule_info['type'] = 'solvent'
        
        self.logger.info(f"处理分子 {molecule_name}，类别: {category}")
        
        # 创建分子子目录
        sanitized_name = molecule_name.replace('+', '').replace('-', '')
        molecule_subdir = os.path.join(molecule_dir, sanitized_name)
        os.makedirs(molecule_subdir, exist_ok=True)
        self.logger.info(f"创建分子目录: {molecule_subdir}")
        
        # 确认分子目录存在
        if not os.path.exists(molecule_subdir):
            error_msg = f"无法创建分子目录: {molecule_subdir}"
            self.logger.error(error_msg)
            raise IOError(error_msg)
        
        # 检查是否使用分子库
        if hasattr(self, 'use_molecule_library') and self.use_molecule_library and self.molecule_library:
            self.logger.info(f"使用分子库获取分子 {molecule_name}")
            
            # 如果是溶剂且提供了SMILE，先尝试通过SMILE查找
            existing_name = None
            if category == "solvents" and molecule_smile:
                existing_name = self.molecule_library.find_by_smile(molecule_smile)
                if existing_name:
                    self.logger.info(f"分子库中找到了匹配SMILE的分子: {existing_name}")
                    
                    # 如果找到的分子名称不同于当前分子名称，使用当前名称获取
                    if existing_name != molecule_name:
                        self.logger.info(f"使用匹配的分子 {existing_name} 代替 {molecule_name}")
            
            # 直接从分子库复制文件（不使用符号链接）
            try:
                # 从分子库获取文件 - 使用复制模式而非符号链接
                source_name = existing_name if existing_name else molecule_name
                self.logger.info(f"从分子库获取分子文件: {source_name} -> {molecule_subdir}")
                
                # 获取源文件路径
                mol_dir = self.molecule_library._get_molecule_path(source_name)
                self.logger.info(f"分子库中的源目录: {mol_dir}")
                
                if not os.path.exists(mol_dir):
                    self.logger.error(f"分子库中的源目录不存在: {mol_dir}")
                    # 继续尝试其他方法
                else:
                    # 手动复制文件
                    result_files = {}
                    
                    # 确定要复制哪些文件扩展名
                    extensions = ["pdb", "lt"] if category in ["cations", "anions"] else ["pdb", "lt", "chg"]
                    
                    for ext in extensions:
                        src_file = os.path.join(mol_dir, f"{source_name}.{ext}")
                        dst_file = os.path.join(molecule_subdir, f"{sanitized_name}.{ext}")
                        
                        self.logger.info(f"检查源文件是否存在: {src_file}")
                        if not os.path.exists(src_file):
                            self.logger.warning(f"源文件不存在: {src_file}")
                            continue
                            
                        # 验证源文件可读性
                        try:
                            with open(src_file, 'rb') as test_file:
                                # 读取一些数据以确保文件可读
                                data = test_file.read(100)
                                self.logger.info(f"源文件 {src_file} 可读，大小: {len(data)} 字节")
                        except Exception as e:
                            self.logger.error(f"源文件不可读: {src_file}, 错误: {str(e)}")
                            continue
                        
                        # 复制文件
                        try:
                            shutil.copy2(src_file, dst_file)
                            self.logger.info(f"已复制文件: {src_file} -> {dst_file}")
                            
                            # 验证目标文件存在并可读
                            if os.path.exists(dst_file):
                                with open(dst_file, 'rb') as test_file:
                                    # 读取一些数据以确保文件可读
                                    data = test_file.read(100)
                                    self.logger.info(f"目标文件 {dst_file} 已创建并可读，大小: {len(data)} 字节")
                                result_files[ext] = dst_file
                            else:
                                self.logger.error(f"目标文件创建失败: {dst_file}")
                        except Exception as e:
                            self.logger.error(f"复制文件时出错: {src_file} -> {dst_file}, 错误: {str(e)}")
                    
                    # 如果成功获取了文件，返回结果
                    if result_files:
                        self.logger.info(f"从分子库成功复制 {source_name} 文件: {result_files}")
                        return result_files
                    else:
                        self.logger.warning(f"无法从分子库获取有效的 {source_name} 文件")
            except Exception as e:
                self.logger.error(f"从分子库获取 {molecule_name} 文件时出错: {str(e)}")
                # 继续尝试其他方法
        
        # 如果这是阳离子或阴离子，尝试从初始盐目录获取
        if category in ["cations", "anions"]:
            if category == "cations":
                salt_dir = self.initial_salts_cations_dir
            else:
                salt_dir = self.initial_salts_anions_dir
            
            if salt_dir and os.path.exists(salt_dir):
                self.logger.info(f"尝试从初始离子目录获取 {molecule_name} 文件: {salt_dir}")
                try:
                    # 从初始离子目录复制文件
                    files = self._copy_existing_molecule_files(molecule_name, salt_dir, molecule_subdir)
                    if files:
                        self.logger.info(f"成功从初始离子目录复制 {molecule_name} 文件: {files}")
                        return files
                    else:
                        self.logger.warning(f"初始离子目录中没有找到 {molecule_name} 文件")
                except Exception as e:
                    self.logger.error(f"从初始离子目录复制 {molecule_name} 文件时出错: {str(e)}")
            else:
                self.logger.warning(f"初始离子目录不存在: {salt_dir}")
        
        # 如果这是溶剂，尝试从初始溶剂目录获取
        if category == "solvents":
            if self.initial_solvent_dir and os.path.exists(self.initial_solvent_dir):
                self.logger.info(f"尝试从初始溶剂目录获取 {molecule_name} 文件: {self.initial_solvent_dir}")
                try:
                    # 从初始溶剂目录复制文件
                    files = self._copy_existing_molecule_files(molecule_name, self.initial_solvent_dir, molecule_subdir)
                    if files:
                        self.logger.info(f"成功从初始溶剂目录复制 {molecule_name} 文件: {files}")
                        return files
                    else:
                        self.logger.warning(f"初始溶剂目录中没有找到 {molecule_name} 文件")
                except Exception as e:
                    self.logger.error(f"从初始溶剂目录复制 {molecule_name} 文件时出错: {str(e)}")
            else:
                self.logger.warning(f"初始溶剂目录不存在: {self.initial_solvent_dir}")
        
            # 如果提供了SMILE，尝试使用LigParGen生成
            if molecule_smile:
                self.logger.info(f"尝试使用LigParGen生成 {molecule_name} 文件")
                try:
                    files = self._generate_molecule_files_with_ligpargen(molecule_info, molecule_dir)
                    if files:
                        self.logger.info(f"成功使用LigParGen生成 {molecule_name} 文件: {files}")
                        
                        # 如果生成了文件，将其添加到分子库
                        if files and hasattr(self, 'use_molecule_library') and self.use_molecule_library and self.molecule_library:
                            smile = molecule_info.get('smile')
                            if smile:
                                try:
                                    self.molecule_library.add_molecule(molecule_name, files, "solvents", smile)
                                    self.logger.info(f"已将溶剂 {molecule_name} 添加到分子库")
                                except Exception as e:
                                    self.logger.error(f"将溶剂 {molecule_name} 添加到分子库时出错: {str(e)}")
                        
                        return files
                    else:
                        self.logger.warning(f"使用LigParGen生成 {molecule_name} 文件失败")
                except Exception as e:
                    self.logger.error(f"使用LigParGen生成 {molecule_name} 文件时出错: {str(e)}")
            else:
                self.logger.warning(f"溶剂 {molecule_name} 没有提供SMILE字符串，无法使用LigParGen生成")
        
        # 所有尝试都失败
        self.logger.error(f"无法生成或获取 {molecule_name} 文件")
        return {}

    def _copy_existing_molecule_files(self, molecule_name, src_dir, dest_dir):
        """从初始目录复制现有分子文件
        
        Args:
            molecule_name: 分子名称
            src_dir: 源目录
            dest_dir: 目标目录
            
        Returns:
            Dict: 包含复制文件路径的字典，如果失败则返回空字典
        """
        self.logger.info(f"从 {src_dir} 复制 {molecule_name} 文件到 {dest_dir}")
        
        # 检查源目录存在
        if not os.path.exists(src_dir):
            self.logger.error(f"源目录不存在: {src_dir}")
            return {}
            
        # 确保目标目录存在
        os.makedirs(dest_dir, exist_ok=True)
        if not os.path.exists(dest_dir):
            self.logger.error(f"无法创建目标目录: {dest_dir}")
            return {}
        
        # 处理不同扩展名的文件
        result_files = {}
        
        # 离子或溶剂可能有的文件扩展名
        extensions = ["pdb", "lt", "chg", "lmp"]
        
        for ext in extensions:
            # 查找源文件 - 考虑大小写变化
            src_files = []
            
            # 标准名称
            standard_src_file = os.path.join(src_dir, f"{molecule_name}.{ext}")
            if os.path.exists(standard_src_file):
                src_files.append(standard_src_file)
                
            # 小写扩展名
            lower_src_file = os.path.join(src_dir, f"{molecule_name}.{ext.lower()}")
            if os.path.exists(lower_src_file) and lower_src_file != standard_src_file:
                src_files.append(lower_src_file)
                
            # 大写扩展名
            upper_src_file = os.path.join(src_dir, f"{molecule_name}.{ext.upper()}")
            if os.path.exists(upper_src_file) and upper_src_file != standard_src_file:
                src_files.append(upper_src_file)
                
            # 如果找到了文件
            if src_files:
                src_file = src_files[0]  # 使用第一个找到的文件
                dst_file = os.path.join(dest_dir, f"{molecule_name}.{ext}")
                
                # 检查源文件是否可读
                try:
                    with open(src_file, 'rb') as f:
                        data = f.read(100)  # 读取一些数据以确保文件可读
                        self.logger.info(f"源文件 {src_file} 可读，大小: {len(data)} 字节")
                except Exception as e:
                    self.logger.error(f"源文件不可读: {src_file}, 错误: {str(e)}")
                    continue
                
                # 复制文件
                try:
                    shutil.copy2(src_file, dst_file)
                    self.logger.info(f"已复制文件: {src_file} -> {dst_file}")
                    
                    # 检查目标文件是否创建成功
                    if os.path.exists(dst_file):
                        # 验证目标文件可读性
                        try:
                            with open(dst_file, 'rb') as f:
                                data = f.read(100)  # 读取一些数据以确保文件可读
                                self.logger.info(f"目标文件 {dst_file} 已创建并可读，大小: {len(data)} 字节")
                            result_files[ext] = dst_file
                        except Exception as e:
                            self.logger.error(f"目标文件不可读: {dst_file}, 错误: {str(e)}")
                    else:
                        self.logger.error(f"目标文件未创建成功: {dst_file}")
                except Exception as e:
                    self.logger.error(f"复制文件时出错: {src_file} -> {dst_file}, 错误: {str(e)}")
        
        return result_files

    def _generate_molecule_files_with_ligpargen(self, molecule_info, molecule_dir):
        """使用LigParGen生成溶剂分子文件
        
        Args:
            molecule_info: 分子信息字典，包含name和smile
            molecule_dir: 输出分子目录
            
        Returns:
            Dict: 包含生成文件路径的字典，如果失败则返回空字典
        """
        # 打印版本标识，帮助验证使用了哪个版本的代码
        self.logger.info("=== 使用修改后的LigParGen处理函数版本3.0 ===")
        
        molecule_name = molecule_info.get('name', '')
        molecule_smile = molecule_info.get('smile', '')
        
        if not molecule_name or not molecule_smile:
            self.logger.warning(f"分子名称或SMILE字符串为空，无法使用LigParGen生成: {molecule_name}")
            return {}
        
        # 首先通过SMILE字符串检查分子库中是否已有该分子
        if hasattr(self, 'use_molecule_library') and self.use_molecule_library and self.molecule_library:
            existing_name = self.molecule_library.find_by_smile(molecule_smile)
            if existing_name:
                self.logger.info(f"在分子库中找到匹配SMILE字符串的分子: {existing_name}，使用该分子代替 {molecule_name}")
                
                # 创建分子子目录
                sanitized_name = molecule_name.replace('+', '').replace('-', '')
                molecule_subdir = os.path.join(molecule_dir, sanitized_name)
                os.makedirs(molecule_subdir, exist_ok=True)
                
                # 从分子库获取文件
                result_files = self.molecule_library.get_molecule(existing_name, molecule_subdir, use_symlink=True)
                
                if result_files:
                    self.logger.info(f"从分子库成功获取匹配的分子 {existing_name} 文件: {result_files}")
                    
                    # 重命名文件以匹配当前分子名称
                    renamed_files = {}
                    for ext, file_path in result_files.items():
                        new_file_path = os.path.join(molecule_subdir, f"{sanitized_name}.{ext}")
                        if os.path.exists(new_file_path):
                            os.remove(new_file_path)
                        
                        # 复制而不是重命名，以保持原始文件不变
                        shutil.copy2(file_path, new_file_path)
                        renamed_files[ext] = new_file_path
                        self.logger.info(f"重命名文件: {file_path} -> {new_file_path}")
                    
                    return renamed_files
        
        self.logger.info(f"使用LigParGen生成溶剂分子 {molecule_name}")
        
        # 创建分子子目录
        sanitized_name = molecule_name.replace('+', '').replace('-', '')
        molecule_subdir = os.path.join(molecule_dir, sanitized_name)
        os.makedirs(molecule_subdir, exist_ok=True)
        
        # 尝试从已有溶剂目录复制
        initial_solvent_dir = self.initial_solvent_path
        if os.path.exists(initial_solvent_dir):
            src_files = self._copy_existing_molecule_files(molecule_name, initial_solvent_dir, molecule_dir)
            if src_files:
                self.logger.info(f"从初始溶剂目录成功复制 {molecule_name} 的文件")
                
                # 将成功生成的文件添加到分子库
                if hasattr(self, 'use_molecule_library') and self.use_molecule_library and self.molecule_library:
                    self.molecule_library.add_molecule(molecule_name, src_files, "solvents", molecule_smile)
                    self.logger.info(f"已将溶剂 {molecule_name} 添加到分子库，SMILE: {molecule_smile}")
                
                return src_files
            
        # 使用ligpargen路径
        ligpargen_path = self.ligpargen_path
        if not ligpargen_path or not os.path.exists(ligpargen_path):
            self.logger.error(f"LigParGen路径不存在: {ligpargen_path}")
            return {}
        
        try:
            # 切换到输出目录
            original_dir = os.getcwd()
            os.chdir(molecule_dir)
            
            # 准备输出文件路径
            pdb_file = os.path.join(molecule_subdir, f"{sanitized_name}.pdb")
            lmp_file = os.path.join(molecule_subdir, f"{sanitized_name}.lmp")
            lt_file = os.path.join(molecule_subdir, f"{sanitized_name}.lt")
            
            # 构建LigParGen命令
            cmd_str = f"{ligpargen_path} -s '{molecule_smile}' -n {sanitized_name} -r MOL -c 0 -o 0 -cgen CM1A"
            
            self.logger.info(f"执行LigParGen命令: {cmd_str}")
            
            # 执行命令
            process = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            
            # 检查命令执行是否成功
            if process.returncode != 0:
                error_msg = f"LigParGen命令执行失败: {stderr}"
                self.logger.error(error_msg)
                if stdout:
                    self.logger.debug(f"标准输出: {stdout}")
                raise RuntimeError(error_msg)
            
            # 检查输出文件是否生成
            result_files = {}
            
            # 列出当前目录的文件
            current_files = os.listdir(os.getcwd())
            self.logger.info(f"当前目录文件列表 (修改版本2.0): {current_files}")
            
            # 检查所有可能的PDB文件名格式
            self.logger.info("开始搜索PDB文件（优先查找.charmm.pdb）...")
            
            # 特定查找.charmm.pdb
            charmm_pdb_file = os.path.join(os.getcwd(), f"{sanitized_name}.charmm.pdb")
            self.logger.info(f"检查charmm.pdb文件: {charmm_pdb_file}")
            
            # 查找任何.pdb文件
            all_pdb_files = [f for f in current_files if f.endswith('.pdb')]
            self.logger.info(f"所有可能的PDB文件: {all_pdb_files}")
            
            # 优先使用.charmm.pdb
            if os.path.exists(charmm_pdb_file):
                self.logger.info(f"找到charmm.pdb文件: {charmm_pdb_file}")
                # 复制到标准位置
                shutil.copy2(charmm_pdb_file, pdb_file)
                result_files['pdb'] = pdb_file
                self.logger.info(f"已复制 {charmm_pdb_file} 到 {pdb_file}")
            elif os.path.exists(pdb_file):
                result_files['pdb'] = pdb_file
                self.logger.info(f"已找到标准PDB文件: {pdb_file}")
            elif all_pdb_files:
                # 使用找到的第一个.pdb文件
                first_pdb = os.path.join(os.getcwd(), all_pdb_files[0])
                self.logger.info(f"使用发现的PDB文件: {first_pdb}")
                shutil.copy2(first_pdb, pdb_file)
                result_files['pdb'] = pdb_file
                self.logger.info(f"已复制 {first_pdb} 到 {pdb_file}")
            else:
                self.logger.warning(f"未找到任何PDB文件")
            
            # 检查所有可能的LMP文件名格式
            self.logger.info("开始搜索LMP文件（优先查找.lammps.lmp）...")
            
            # 特定查找.lammps.lmp
            lammps_lmp_file = os.path.join(os.getcwd(), f"{sanitized_name}.lammps.lmp")
            self.logger.info(f"检查lammps.lmp文件: {lammps_lmp_file}")
            
            # 查找任何.lmp文件
            all_lmp_files = [f for f in current_files if f.endswith('.lmp')]
            self.logger.info(f"所有可能的LMP文件: {all_lmp_files}")
            
            # 优先使用.lammps.lmp
            if os.path.exists(lammps_lmp_file):
                self.logger.info(f"找到lammps.lmp文件: {lammps_lmp_file}")
                # 复制到标准位置
                shutil.copy2(lammps_lmp_file, lmp_file)
                result_files['lmp'] = lmp_file
                self.logger.info(f"已复制 {lammps_lmp_file} 到 {lmp_file}")
            elif os.path.exists(lmp_file):
                result_files['lmp'] = lmp_file
                self.logger.info(f"已找到标准LMP文件: {lmp_file}")
            elif all_lmp_files:
                # 使用找到的第一个.lmp文件
                first_lmp = os.path.join(os.getcwd(), all_lmp_files[0])
                self.logger.info(f"使用发现的LMP文件: {first_lmp}")
                shutil.copy2(first_lmp, lmp_file)
                result_files['lmp'] = lmp_file
                self.logger.info(f"已复制 {first_lmp} 到 {lmp_file}")
            else:
                self.logger.warning(f"未找到任何LMP文件")
            
            # 检查是否需要生成电荷文件
            if 'pdb' in result_files and self.use_resp_charges and self.resp_path:
                charge_file = os.path.join(molecule_subdir, f"{sanitized_name}.chg")
                # 如果还没有电荷文件，尝试使用RESP生成
                if not os.path.exists(charge_file) and os.path.exists(self.resp_path):
                    self.logger.info(f"尝试使用RESP2为分子{sanitized_name}生成电荷文件")
                    try:
                        # 调用RESP生成电荷
                        self._generate_charges_with_resp(sanitized_name, result_files['pdb'], self.resp_path)
                        
                        # 检查是否生成了电荷文件
                        resp_charge_file = os.path.join(self.resp_path, f"{sanitized_name}.charge")
                        if os.path.exists(resp_charge_file):
                            self.logger.info(f"RESP生成了电荷文件: {resp_charge_file}")
                            # 复制并重命名为.chg文件
                            shutil.copy2(resp_charge_file, charge_file)
                            self.logger.info(f"已复制电荷文件: {resp_charge_file} -> {charge_file}")
                            
                            # 如果有LMP文件，更新其中的电荷
                            if 'lmp' in result_files:
                                try:
                                    self.logger.info(f"使用RESP电荷更新LMP文件中的电荷")
                                    replace_charges_in_lmp(result_files['lmp'], charge_file)
                                    self.logger.info(f"成功更新LMP文件中的电荷")
                                except Exception as e:
                                    self.logger.error(f"更新LMP文件电荷时出错: {str(e)}")
                    except Exception as e:
                        self.logger.error(f"使用RESP生成电荷时出错: {str(e)}")
            
            # 检查或生成LT文件
            if os.path.exists(lt_file):
                result_files['lt'] = lt_file
                self.logger.info(f"已找到LT文件: {lt_file}")
            elif 'lmp' in result_files:
                # 如果LMP文件存在但LT文件不存在，尝试从LMP生成LT
                try:
                    self.logger.info(f"尝试从LMP文件生成LT文件")
                    lt_file = self._generate_lt_from_lmp(result_files['lmp'], lt_file, sanitized_name)
                    if os.path.exists(lt_file):
                        result_files['lt'] = lt_file
                        self.logger.info(f"成功从LMP生成LT文件: {lt_file}")
                except Exception as e:
                    self.logger.error(f"从LMP生成LT文件失败: {e}")
            
            # 总结找到的文件
            self.logger.info(f"文件生成结果: PDB: {'是' if 'pdb' in result_files else '否'}, LMP: {'是' if 'lmp' in result_files else '否'}, LT: {'是' if 'lt' in result_files else '否'}")
            
            return result_files
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"LigParGen执行失败: {e}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            
            # 记录错误日志
            error_log = os.path.join(molecule_dir, "ligpargen_error.log")
            with open(error_log, 'w') as f:
                f.write(f"Command: {cmd_str}\n")
                f.write(f"Error: {str(e)}\n")
                if e.stdout:
                    f.write("\n--- STDOUT ---\n")
                    f.write(e.stdout)
                if e.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(e.stderr)
            
            raise RuntimeError(f"LigParGen执行失败，详情查看: {error_log}")
        except Exception as e:
            self.logger.error(f"使用LigParGen生成分子文件时出错: {e}")
            return {}
            
    def _generate_lt_from_lmp(self, lmp_file, lt_file, molecule_name):
        """使用ltemplify从LMP文件生成LT文件
        
        Args:
            lmp_file: LMP文件路径
            lt_file: 输出LT文件路径
            molecule_name: 分子名称
            
        Returns:
            str: 生成的LT文件路径
            
        Raises:
            FileNotFoundError: 当LMP文件不存在或ltemplify命令不可用时
            RuntimeError: 当ltemplify执行失败时
        """
        self.logger.info(f"使用ltemplify生成LT文件: {lt_file}")
        
        # 检查LMP文件是否存在
        if not os.path.exists(lmp_file):
            error_msg = f"LMP文件不存在: {lmp_file}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
            
        # 查找ltemplify命令
        ltemplify_cmd = self.ltemplify_path
        
        # 确保ltemplify命令是一个存在的文件
        if not os.path.exists(ltemplify_cmd):
            # 尝试其他可能的路径
            possible_paths = [
                "/public/software/moltemplate_2025-2-02/moltemplate/ltemplify.py",
                "/public/software/moltemplate/moltemplate/ltemplify.py",
                os.path.expanduser("~/bin/ltemplify.py")
            ]
            
            for path in possible_paths:
                if os.path.exists(path):
                    ltemplify_cmd = path
                    self.logger.info(f"找到ltemplify文件: {ltemplify_cmd}")
                    break
                    
            if not os.path.exists(ltemplify_cmd):
                error_msg = f"找不到ltemplify.py文件，尝试过的路径: {[ltemplify_cmd] + possible_paths}"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
                
        # 确保输出目录存在
        os.makedirs(os.path.dirname(lt_file), exist_ok=True)
        
        # 执行ltemplify命令
        python_cmd = "python"
        cmd = f"{python_cmd} {ltemplify_cmd} -name {molecule_name} {lmp_file} > {lt_file}"
        
        self.logger.info(f"执行ltemplify命令: {cmd}")
        
        try:
            # 使用subprocess执行命令
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            
            # 检查命令执行是否成功
            if process.returncode != 0:
                error_msg = f"ltemplify命令执行失败: {stderr}"
                self.logger.error(error_msg)
                if stdout:
                    self.logger.debug(f"标准输出: {stdout}")
                raise RuntimeError(error_msg)
                
            # 检查生成的LT文件
            if os.path.exists(lt_file):
                # 检查文件是否有内容
                with open(lt_file, 'r') as f:
                    content = f.read().strip()
                    
                if not content:
                    error_msg = f"ltemplify生成的LT文件内容为空: {lt_file}"
                    self.logger.error(error_msg)
                    raise RuntimeError(error_msg)
                    
                self.logger.info(f"成功生成LT文件: {lt_file}")
                return lt_file
            else:
                error_msg = f"ltemplify未能生成LT文件: {lt_file}"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
                
        except Exception as e:
            self.logger.error(f"执行ltemplify时出错: {e}")
            raise

    def _generate_lammps_input_with_moltemplate(self, system_lt_file, output_dir):
        """使用moltemplate生成LAMMPS输入文件
        
        Args:
            system_lt_file: 系统LT文件路径
            output_dir: 输出目录
            
        Returns:
            dict: 包含生成的文件路径的字典
        """
        self.logger.info(f"使用moltemplate生成LAMMPS输入文件: {system_lt_file}")
        
        # 检查系统LT文件是否存在
        if not os.path.exists(system_lt_file):
            error_msg = f"系统LT文件不存在: {system_lt_file}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        # 检查moltemplate路径是否存在
        moltemplate_path = self.moltemplate_path
        if not moltemplate_path or not os.path.exists(moltemplate_path):
            error_msg = f"moltemplate路径不存在: {moltemplate_path}"
            self.logger.error(error_msg)
            raise FileNotFoundError(error_msg)
            
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # 首先收集所有所需的LT文件
            lt_files_dir = os.path.dirname(system_lt_file)
            self.logger.info(f"从 {lt_files_dir} 收集LT文件到 {output_dir}")
            
            # 分析system.lt文件找出所有导入的LT文件
            imported_files = []
            with open(system_lt_file, 'r') as f:
                for line in f:
                    if line.strip().startswith('import'):
                        # 提取导入文件名
                        match = line.strip().split('"')[1]
                        imported_files.append(match)
            self.logger.info(f"系统LT文件导入的文件: {imported_files}")
            
            # 1. 复制所有分子目录下的.lt文件到输出目录
            source_dirs = []
            for name, files in self.molecule_files.items():
                if 'lt' in files:
                    lt_dir = os.path.dirname(files['lt'])
                    if lt_dir not in source_dirs:
                        source_dirs.append(lt_dir)
            
            # 从各个源目录复制lt文件
            for src_dir in source_dirs:
                for lt_file in glob.glob(os.path.join(src_dir, "*.lt")):
                    lt_name = os.path.basename(lt_file)
                    dest_lt = os.path.join(output_dir, lt_name)
                    shutil.copy2(lt_file, dest_lt)
                    self.logger.info(f"已复制LT文件: {lt_file} -> {dest_lt}")
            
            # 2. 复制系统LT文件到输出目录
            dest_system_lt = os.path.join(output_dir, "system.lt")
            if system_lt_file != dest_system_lt:
                shutil.copy2(system_lt_file, dest_system_lt)
                self.logger.info(f"已复制系统LT文件: {system_lt_file} -> {dest_system_lt}")
            
            # 3. 验证所有导入的文件是否存在
            missing_files = []
            for imp_file in imported_files:
                if not os.path.exists(os.path.join(output_dir, imp_file)):
                    missing_files.append(imp_file)
            
            if missing_files:
                self.logger.warning(f"以下导入的LT文件在输出目录中不存在: {missing_files}")
                
                # 尝试查找这些文件并复制
                for missing in missing_files:
                    # 在整个工作目录中搜索
                    found_files = []
                    for root, _, files in os.walk(self.work_dir):
                        if missing in files:
                            found_files.append(os.path.join(root, missing))
                    
                    if found_files:
                        # 使用找到的第一个文件
                        shutil.copy2(found_files[0], os.path.join(output_dir, missing))
                        self.logger.info(f"找到并复制了缺失的文件: {found_files[0]} -> {os.path.join(output_dir, missing)}")
                    else:
                        self.logger.error(f"无法找到缺失的文件: {missing}")
            
            # 保存当前工作目录
            original_dir = os.getcwd()
            
            # 切换到输出目录
            os.chdir(output_dir)
            self.logger.info(f"切换到目录: {output_dir}")
            
            # 再次验证所有输出目录下的文件
            existing_files = os.listdir(output_dir)
            self.logger.info(f"执行moltemplate前输出目录中的文件: {existing_files}")
            
            # 构建命令 - 使用绝对路径指定moltemplate，但处理相对于当前目录的system.lt
            cmd = f"{moltemplate_path} -pdb system.pdb system.lt"
            self.logger.info(f"执行moltemplate命令: {cmd}")
            
            # 执行命令
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            
            # 记录完整的输出
            self.logger.info(f"Moltemplate stdout: {stdout}")
            if stderr:
                self.logger.warning(f"Moltemplate stderr: {stderr}")
            
            # 检查命令执行是否成功
            if process.returncode != 0:
                error_msg = f"moltemplate命令执行失败: {stderr}"
                self.logger.error(error_msg)
                # 写入错误日志文件
                with open(os.path.join(output_dir, "moltemplate_error.log"), 'w') as f:
                    f.write(f"命令: {cmd}\n")
                    f.write(f"返回码: {process.returncode}\n")
                    f.write(f"标准输出:\n{stdout}\n")
                    f.write(f"错误输出:\n{stderr}\n")
                raise RuntimeError(error_msg)
            
            self.logger.info(f"moltemplate命令执行成功")
                
            # 收集生成的文件
            generated_files = {
                'data_file': os.path.join(output_dir, "system.data"),
                'in_file': os.path.join(output_dir, "system.in"),
                'in_settings_file': os.path.join(output_dir, "system.in.settings"),
                'in_init_file': os.path.join(output_dir, "system.in.init"),
                'pdb_file': os.path.join(output_dir, "system.pdb")
            }
            
            # 检查文件是否存在
            for file_key, file_path in list(generated_files.items()):
                if not os.path.exists(file_path):
                    self.logger.warning(f"{file_key}文件不存在: {file_path}")
                    # 从结果中移除不存在的文件
                    del generated_files[file_key]
                else:
                    self.logger.info(f"生成的{file_key}文件: {file_path}")
            
            # 恢复原始目录
            if 'original_dir' in locals():
                os.chdir(original_dir)
                self.logger.info(f"恢复原始目录: {original_dir}")
                
            # 生成LAMMPS输入脚本
            self._generate_lammps_input_script(system_lt_file, output_dir)
            
            return generated_files
            
        except Exception as e:
            self.logger.error(f"执行moltemplate时发生异常: {str(e)}", exc_info=True)
            raise RuntimeError(f"执行moltemplate时发生异常: {str(e)}")
            
        finally:
            # 恢复原始工作目录
            if 'original_dir' in locals():
                os.chdir(original_dir)
                self.logger.info(f"恢复原始目录: {original_dir}")

    def _generate_lammps_input_script(self, system_lt_file, output_dir):
        """生成LAMMPS输入脚本
        
        Args:
            system_lt_file: 系统LT文件路径
            output_dir: 输出目录
            
        Returns:
            str: 生成的LAMMPS输入文件路径
        """
        self.logger.info(f"生成LAMMPS输入脚本")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 从formulation_data['parameters']获取参数，使用默认值
        parameters = self.formulation_data.get('parameters', {})
        temp = parameters.get('temperature', 298.15)
        pressure = parameters.get('pressure', 1.0)
        cutoff = parameters.get('cutoff', 10.0)
        time_step = parameters.get('time_step', 1.0)
        box_size = float(parameters.get('box_size', 50.0))
        equilibration_steps = parameters.get('equilibration_steps', 10000)
        production_steps = parameters.get('production_steps', 100000)
        
        # 设置输出频率
        thermo_freq = 1000
        trj_freq_npt = 2000
        trj_freq_nvt = 10000
        
        # 确定系统名称
        system_name = os.path.basename(system_lt_file).replace('.lt', '')
        
        # 获取系统成分信息
        cations = self.config.get('cations', [])
        anions = self.config.get('anions', [])
        solvents = self.config.get('solvents', [])
        
        # 创建元素列表
        # 这里需要根据实际情况调整，或者从LAMMPS数据文件中读取
        elements = []
        for cation in cations:
            if cation.get('name') == 'Li':
                elements.append('Li')
        for anion in anions:
            if anion.get('name') == 'PF6':
                elements.append('P')
                elements.extend(['F'] * 6)
        for solvent in solvents:
            if solvent.get('name') == 'EC':
                # EC分子大致包含以下元素
                elements.extend(['C'] * 3)
                elements.extend(['O'] * 2)
                elements.extend(['H'] * 4)
                
        element_list = " ".join(elements) if elements else "Li P F F F F F F C C C O O H H H H"
        
        # 创建RDF对的信息
        rdf_pairs = ""
        rdf_titles = ""
        
        # Li-F RDF
        if any(cat.get('name') == 'Li' for cat in cations) and any(an.get('name') == 'PF6' for an in anions):
            rdf_pairs += "1 3 1 4 1 5 1 6 1 7 1 8 "
            rdf_titles += "rdf_Li_F "
        
        # Li-O RDF
        if any(cat.get('name') == 'Li' for cat in cations) and any(sol.get('name') == 'EC' for sol in solvents):
            rdf_pairs += "1 11 1 12 "
            rdf_titles += "rdf_Li_O "
            
        # 生成LAMMPS输入文件内容
        in_file_content = f"""# LAMMPS输入脚本 - 生成于 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# ----------------- Variable Section -----------------
variable infile string {system_name}
variable outname string {system_name}

include {system_name}.in.init
read_data {system_name}.data
include {system_name}.in.settings

# 元素列表
variable element_list index "{element_list}"
variable rdf_pair string "{rdf_pairs}"

# 模拟参数
variable Nsteps_NPT equal {equilibration_steps}
variable Nsteps_NVT equal {production_steps}
variable Freq_trj_npt equal {trj_freq_npt}
variable Freq_trj_nvt equal {trj_freq_nvt}
variable thermo_freq equal {thermo_freq}
variable timestep equal {time_step}
variable Temp_NPT equal {temp}
variable Temp_NVT equal {temp}

# 创建分子组
"""
        
        # 添加分子组信息
        for i, cat in enumerate(cations, start=1):
            if 'name' in cat and cat.get('number', 0) > 0:
                cat_name = cat['name'].replace('+', '').replace('-', '')
                in_file_content += f'group Cation{i} type {i}\n'
                
        for i, anion in enumerate(anions, start=1):
            if 'name' in anion and anion.get('number', 0) > 0:
                anion_name = anion['name'].replace('+', '').replace('-', '')
                in_file_content += f'group Anion{i} type {i+len(cations)}\n'
                
        # 主模拟部分
        in_file_content += f"""
# 输出设置
thermo_style custom step cpu cpuremain temp density lx ly lz etotal ke pe evdwl ecoul elong ebond eangle edihed eimp
thermo ${{thermo_freq}}
timestep ${{timestep}}

# 能量最小化
minimize 1.0e-4 1.0e-6 5000 10000

# 保存最小化后的结构
write_data ${{infile}}_after_minimize.data nocoeff
write_dump all custom ${{infile}}_after_minimize.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id

reset_timestep 0

# NPT平衡阶段
dump trj_npt all custom ${{Freq_trj_npt}} NPT_${{outname}}.lammpstrj id element mol type x y z q
dump_modify trj_npt flush yes element ${{element_list}} sort id

fix fxnpt all npt temp ${{Temp_NPT}} ${{Temp_NPT}} $(100.0*dt) iso {pressure} {pressure} $(1000.0*dt)
run ${{Nsteps_NPT}}
unfix fxnpt

undump trj_npt

# 保存NPT后的结构
write_data ${{infile}}_after_npt.data
write_restart ${{infile}}_restart_after_npt.data
write_dump all custom ${{infile}}_after_npt.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id

reset_timestep 0

# 计算MSD
"""

        # 添加MSD计算
        for i, cation in enumerate(cations, start=1):
            if 'name' in cation and cation.get('number', 0) > 0:
                in_file_content += f'compute Ca{i} Cation{i} msd com yes\n'
                in_file_content += f'fix Ca{i}msd Cation{i} ave/time ${{Freq_trj_nvt}} 1 ${{Freq_trj_nvt}} c_Ca{i}[1] c_Ca{i}[2] c_Ca{i}[3] c_Ca{i}[4] file out_cation{i}_msd.dat title1 "t msd msd msd msd_cation{i}" title2 "fs x y z total"\n'
                
        for i, anion in enumerate(anions, start=1):
            if 'name' in anion and anion.get('number', 0) > 0:
                in_file_content += f'compute An{i} Anion{i} msd com yes\n'
                in_file_content += f'fix An{i}msd Anion{i} ave/time ${{Freq_trj_nvt}} 1 ${{Freq_trj_nvt}} c_An{i}[1] c_An{i}[2] c_An{i}[3] c_An{i}[4] file out_anion{i}_msd.dat title1 "t msd msd msd msd_anion{i}" title2 "fs x y z total"\n'
                
        # 添加RDF计算和NVT模拟
        in_file_content += f"""
# 计算径向分布函数RDF
compute rdfc1 all rdf 100 ${{rdf_pair}}
fix rdff1 all ave/time $(v_Nsteps_NVT/1000) 1000 ${{Nsteps_NVT}} c_rdfc1[*] file out_rdf.dat mode vector title3 "RDF {rdf_titles}"

# NVT生产阶段
dump trj_nvt all custom ${{Freq_trj_nvt}} NVT_${{outname}}.lammpstrj id element mol type x y z q
dump_modify trj_nvt flush yes element ${{element_list}} sort id
dump utrj_nvt all custom ${{Freq_trj_nvt}} NVT_${{outname}}_un.lammpstrj id element mol type xu yu zu ix iy iz q
dump_modify utrj_nvt flush yes element ${{element_list}} sort id

fix fxnvt all nvt temp ${{Temp_NVT}} ${{Temp_NVT}} $(100.0*dt)
run ${{Nsteps_NVT}}
unfix fxnvt

undump trj_nvt
undump utrj_nvt

# 保存最终结构
write_data ${{infile}}_after_nvt.data
write_restart ${{infile}}_restart_after_nvt.data
write_dump all custom ${{infile}}_after_nvt.lammpstrj id element mol type x y z q modify element ${{element_list}} sort id
"""

        # 写入LAMMPS输入文件
        in_file_path = os.path.join(output_dir, f"{system_name}.in")
        with open(in_file_path, 'w') as f:
            f.write(in_file_content)
            
        self.logger.info(f"生成的LAMMPS输入脚本: {in_file_path}")
        
        # 生成SLURM作业提交脚本
        self._generate_job_script(system_name, in_file_path, output_dir)
        
        return in_file_path
        
    def _generate_job_script(self, system_name, in_file_path, output_dir):
        """生成SLURM作业提交脚本
        
        Args:
            system_name: 系统名称，用作作业名
            in_file_path: LAMMPS输入文件路径
            output_dir: 输出目录路径
            
        Returns:
            str: 生成的作业脚本路径
        """
        self.logger.info(f"生成SLURM作业提交脚本")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 获取输入文件的基本名称
        in_file_name = os.path.basename(in_file_path)
        job_name = system_name
        
        # 从配置中获取CPU数量和最大运行时间
        ncpus = self.config.get('cpus', 64)
        time_hours = self.config.get('max_hours', 72)
        
        # 创建作业脚本路径
        job_script_path = os.path.join(output_dir, "job.sh")
        
        # 创建作业脚本内容
        job_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=out.dat
#SBATCH --error=err.dat
#SBATCH --partition=cpu  # 使用 cpu 队列
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={ncpus}  # 设置每个任务的CPU核心数
#SBATCH --time={time_hours}:00:00  # 设置最大运行时间

# 进入作业提交目录
cd $SLURM_SUBMIT_DIR

EXEC=/public/software/lammps/mpich_install/bin/lmp_mpi
mpirun -n {ncpus} $EXEC < {in_file_name} > {system_name}.log
"""
        
        # 写入作业脚本文件
        with open(job_script_path, 'w') as f:
            f.write(job_script_content)
        
        # 设置执行权限
        os.chmod(job_script_path, 0o755)
        
        self.logger.info(f"生成的SLURM作业提交脚本: {job_script_path}")
        
        # 自动提交作业（如果配置允许）
        if self.config.get('auto_submit', False):
            self._submit_job(job_script_path, output_dir)
        
        return job_script_path
        
    def _submit_job(self, job_script_path, output_dir):
        """提交SLURM作业
        
        Args:
            job_script_path: 作业脚本路径
            output_dir: 输出目录路径
            
        Returns:
            str: 作业ID
        """
        self.logger.info(f"提交SLURM作业")
        
        # 保存当前目录
        original_dir = os.getcwd()
        job_id = None
        
        try:
            # 切换到输出目录
            os.chdir(output_dir)
            self.logger.info(f"切换到输出目录: {output_dir}")
            
            # 设置环境变量
            env = os.environ.copy()
            # 确保PATH中包含sbatch命令
            if '/usr/bin' not in env.get('PATH', ''):
                env['PATH'] = '/usr/bin:' + env.get('PATH', '')
            if '/usr/local/bin' not in env.get('PATH', ''):
                env['PATH'] = '/usr/local/bin:' + env.get('PATH', '')
            self.logger.info(f"设置环境变量 PATH: {env['PATH']}")
            
            # 提交作业
            self.logger.info(f"执行命令: sbatch {job_script_path}")
            import subprocess
            try:
                submit_result = subprocess.run(
                    ["sbatch", job_script_path],
                    capture_output=True,
                    text=True,
                    check=True,
                    env=env
                )
                
                # 解析作业ID
                stdout = submit_result.stdout.strip()
                self.logger.info(f"sbatch命令输出: {stdout}")
                
                # 典型输出格式: "Submitted batch job 123456"
                if "Submitted batch job" in stdout:
                    job_id = stdout.split()[-1]
                    self.logger.info(f"作业已提交，作业ID: {job_id}")
                    
                    # 保存作业ID到结果字典
                    self.results['slurm_job_id'] = job_id
                else:
                    self.logger.warning(f"无法从输出中解析作业ID: {stdout}")
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"提交作业失败，返回码: {e.returncode}")
                self.logger.error(f"错误输出: {e.stderr}")
                self.logger.error(f"标准输出: {e.stdout}")
                
                # 尝试直接执行命令行而不使用subprocess.run
                try:
                    self.logger.info("尝试使用os.system提交作业")
                    status = os.system(f"sbatch {job_script_path}")
                    if status == 0:
                        self.logger.info("使用os.system成功提交作业")
                        job_id = "unknown_id"  # 无法获取确切的作业ID
                        
                        # 尝试从输出文件中获取作业ID
                        try:
                            # 检查当前目录下是否有sbatch输出
                            import glob
                            slurm_output_files = glob.glob("slurm-*.out")
                            if slurm_output_files:
                                # 按时间排序，获取最新的
                                latest_file = max(slurm_output_files, key=os.path.getctime)
                                # 提取作业ID
                                job_id = latest_file.split("-")[1].split(".")[0]
                                self.logger.info(f"从输出文件中解析到作业ID: {job_id}")
                                
                                # 保存作业ID到结果字典
                                self.results['slurm_job_id'] = job_id
                        except Exception as e3:
                            self.logger.warning(f"尝试从输出文件解析作业ID失败: {str(e3)}")
                    else:
                        self.logger.error(f"使用os.system提交作业失败，返回码: {status}")
                except Exception as e2:
                    self.logger.error(f"使用os.system提交作业失败: {str(e2)}")
                
            # 保存作业ID到全局变量
            if job_id:
                # 添加到结果字典
                self.results['slurm_job_id'] = job_id
                
                # 尝试更新数据库中的记录
                if hasattr(self, 'update_database') and callable(self.update_database):
                    try:
                        self.update_database('slurm_job_id', job_id)
                        self.logger.info(f"成功更新数据库中的SLURM作业ID: {job_id}")
                    except Exception as e:
                        self.logger.warning(f"更新数据库作业ID失败: {str(e)}")
                
        except Exception as e:
            self.logger.error(f"提交作业失败: {str(e)}")
            # 获取更详细的错误信息
            import traceback
            self.logger.error(f"详细错误信息: {traceback.format_exc()}")
            
        finally:
            # 恢复原始目录
            os.chdir(original_dir)
            self.logger.info(f"恢复原始目录: {original_dir}")
        
        return job_id

    def _build_system_with_packmol(self, molecules_info, output_dir):
        """使用Packmol构建分子系统
        
        Args:
            molecules_info: 分子信息列表或字典，每个元素为带有files和sanitized_name的字典
            output_dir: 输出目录
            
        Returns:
            dict: 包含生成的系统文件路径的字典
        """
        self.logger.info("开始使用Packmol构建分子系统")
        
        # 打印molecules_info类型和内容，帮助调试
        self.logger.info(f"molecules_info类型: {type(molecules_info)}")
        
        # 处理molecules_info参数，确保是正确的数据结构
        if isinstance(molecules_info, str):
            self.logger.warning(f"molecules_info是字符串，需要从formulation_data中获取实际分子信息")
            # 从formulation_data中获取分子信息
            molecules_info = []
            for cat in self.formulation_data.get('cations', []):
                if isinstance(cat, dict):
                    molecules_info.append(cat)
            for an in self.formulation_data.get('anions', []):
                if isinstance(an, dict):
                    molecules_info.append(an)
            for sol in self.formulation_data.get('solvents', []):
                if isinstance(sol, dict):
                    molecules_info.append(sol)
            
            self.logger.info(f"从formulation_data中重新构建了molecules_info，现在包含{len(molecules_info)}个分子")
        
        # 计算盒子尺寸 - 从formulation_data的parameters中获取，确保参数类型统一
        box_size = float(self.formulation_data.get('parameters', {}).get('box_size', 50.0))
        half_box = box_size / 2.0
        
        # 添加更多日志来追踪box_size的值
        self.logger.info(f"【系统构建】在_build_system_with_packmol中使用box_size: {box_size}")
        self.logger.info(f"【参数检查】当前formulation_data['parameters']内容: {self.formulation_data.get('parameters', {})}")
        
        # 计算和记录盒子体积，用于调试
        box_volume_A3 = box_size ** 3  # 体积，单位：埃^3
        box_volume_nm3 = box_volume_A3 / 1000.0  # 体积，单位：纳米^3
        box_volume_L = box_volume_A3 * 1e-27  # 体积，单位：升
        
        self.logger.info(f"盒子尺寸: {box_size} Å")
        self.logger.info(f"盒子体积: {box_volume_A3:.2f} Å³ = {box_volume_nm3:.2f} nm³ = {box_volume_L:.6f} L")
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        # 输出文件路径
        output_pdb = os.path.join(output_dir, "system.pdb")
        
        # 检查是否有分子信息
        if not molecules_info:
            self.logger.error("未提供任何分子信息或无法从config构建分子信息")
            raise ValueError("未提供任何分子信息或无法从config构建分子信息")
        
        # 整理分子数量并进行日志记录
        total_molecules = 0
        cation_molecules = 0
        anion_molecules = 0
        solvent_molecules = 0
        
        # 生成Packmol输入文件
        packmol_input_file = os.path.join(output_dir, "packmol.inp")
        self.logger.info(f"生成Packmol输入文件: {packmol_input_file}")
        
        try:
            with open(packmol_input_file, 'w') as f:
                # 写入Packmol文件头
                f.write("# Packmol input file generated by MolyteWorkflow\n")
                f.write(f"# System name: {self.formulation_data.get('name', 'unnamed')}\n")
                f.write("\n")
                f.write("tolerance 2.0\n")
                f.write("filetype pdb\n")
                f.write(f"output {output_pdb}\n")
                f.write("\n")
                
                # 写入每种分子的信息
                for molecule in molecules_info:
                    name = molecule.get("name", "")
                    mol_type = molecule.get("type", "")
                    files = molecule.get("files", {})
                    count = molecule.get("number", molecule.get("count", 0))
                    
                    # 跳过没有PDB文件的分子
                    if 'pdb' not in files:
                        self.logger.warning(f"分子 {name} 没有PDB文件，跳过")
                        continue
                    
                    # 跳过数量为0的分子
                    if count <= 0:
                        self.logger.warning(f"分子 {name} 的数量为0，跳过")
                        continue
                    
                    # 更新统计信息
                    total_molecules += count
                    if mol_type == 'cation':
                        cation_molecules += count
                    elif mol_type == 'anion':
                        anion_molecules += count
                    else:
                        solvent_molecules += count
                    
                    # 写入分子信息
                    f.write(f"# Molecule: {name}, Type: {mol_type}, Count: {count}\n")
                    f.write(f"structure {files['pdb']}\n")
                    f.write(f"  number {int(count)}\n")  # 将数量转换为整数
                    f.write(f"  inside box 0. 0. 0. {box_size} {box_size} {box_size}\n")
                    f.write("end structure\n\n")
            
            # 记录统计信息
            self.logger.info(f"总分子数: {total_molecules} (阳离子: {cation_molecules}, 阴离子: {anion_molecules}, 溶剂: {solvent_molecules})")
            
            # 运行Packmol
            result = self._execute_packmol(packmol_input_file, output_pdb, output_dir)
            
            if not result or 'error' in result:
                self.logger.error(f"Packmol执行失败: {result.get('error', '未知错误')}")
                raise RuntimeError(f"Packmol执行失败: {result.get('error', '未知错误')}")
            
            # 收集生成的PDB文件
            if not os.path.exists(output_pdb):
                self.logger.error(f"Packmol未能生成输出PDB文件: {output_pdb}")
                raise FileNotFoundError(f"Packmol未能生成输出PDB文件: {output_pdb}")
            
            # 生成系统LT文件
            self.logger.info("生成系统LT文件")
            
            # 收集分子的LT文件
            lt_files = []
            for molecule in molecules_info:
                files = molecule.get("files", {})
                if 'lt' in files:
                    lt_files.append(files['lt'])
            
            # 生成系统LT文件
            system_lt_file = os.path.join(output_dir, "system.lt")
            
            # 构建分子数量字典
            num_molecules = {mol.get("name", ""): mol.get("number", mol.get("count", 0)) for mol in molecules_info}
            
            # 生成系统LT文件
            self._generate_system_lt_file(lt_files, system_lt_file, self.molecule_files, num_molecules)
            
            # 使用moltemplate生成LAMMPS输入文件
            self.logger.info("使用moltemplate生成LAMMPS输入文件")
            lammps_files = self._generate_lammps_input_with_moltemplate(system_lt_file, output_dir)
            
            # 收集并返回所有文件
            return {
                'pdb': output_pdb,
                'lt': system_lt_file,
                'in': lammps_files.get('in_file'),
                'data': lammps_files.get('data_file')
            }
        
        except Exception as e:
            self.logger.error(f"构建系统时出错: {str(e)}", exc_info=True)
            raise RuntimeError(f"构建系统时出错: {str(e)}")

    def _generate_system_lt_file(self, lt_files, system_lt_file, molecule_files, num_molecules):
        """生成系统LT文件
        
        Args:
            lt_files: LT文件列表
            system_lt_file: 目标系统LT文件路径
            molecule_files: 分子文件字典
            num_molecules: 每个分子的数量字典
            
        Returns:
            str: 生成的系统LT文件路径
        """
        self.logger.info(f"生成系统LT文件: {system_lt_file}")
        
        # 创建导入部分
        imports = []
        for lt_file in lt_files:
            # 使用相对路径
            base_name = os.path.basename(lt_file)
            imports.append(f'import "{base_name}"')
        
        # 获取box_size - 确保与_build_system_with_packmol方法使用一致的值
        box_size = float(self.formulation_data.get('parameters', {}).get('box_size', 50.0))
        
        # 添加更多日志来追踪box_size的值
        self.logger.info(f"【LT文件生成】在_generate_system_lt_file中使用box_size: {box_size}")
        self.logger.info(f"【参数检查】当前formulation_data['parameters']内容: {self.formulation_data.get('parameters', {})}")
        
        # 创建分子实例
        molecule_instances = []
        
        # 标记是否有阳离子和阴离子
        has_cation = False
        has_anion = False
        
        # 处理所有分子
        for mol_name, count in num_molecules.items():
            if count <= 0:
                continue
            
            # 处理分子名称，移除+或-符号
            clean_mol_name = mol_name.replace('+', '').replace('-', '')
            
            # 检查是否为阳离子或阴离子
            is_cation = False
            is_anion = False
            
            for cat in self.config.get("cations", []):
                if cat.get("name") == mol_name:
                    is_cation = True
                    has_cation = True
                    break
            
            for an in self.config.get("anions", []):
                if an.get("name") == mol_name:
                    is_anion = True
                    has_anion = True
                    break
            
            # 添加分子实例，使用数组形式，确保数量是整数
            molecule_instances.append(f"{clean_mol_name}s = new {clean_mol_name}[{int(count)}]")
        
        # 写入系统LT文件
        with open(system_lt_file, 'w') as f:
            # 导入分子类型定义
            for imp in imports:
                f.write(f"{imp}\n")
            
            f.write("\n\n\n")
            
            # 添加注释说明Moltemplate用法
            f.write("# In this example we are using PACKMOL so we don't have to specify the position\n")
            f.write("# of every molecule.  Instead we list the number of molecules of each type\n")
            f.write("# that we need.  However you must be careful to list the lipids and water in\n")
            f.write("# the same order they appear in the \"mix_lipids+water.inp\" (PACKMOL input) file\n\n")
            
            # 添加分子实例定义
            for instance in molecule_instances:
                f.write(f"{instance}\n\n")
            
            f.write("\n\n")
            
            # 添加模拟盒子大小定义
            f.write("# How big is the simulation box?\n\n")
            f.write("write_once(\"Data Boundary\") {\n")
            f.write(f"0 {box_size} xlo xhi\n")
            f.write(f"0 {box_size} ylo yhi\n")
            f.write(f"0 {box_size} zlo zhi\n")
            f.write("}\n")
        
        self.logger.info(f"成功生成系统LT文件: {system_lt_file}")
        return system_lt_file

    def _generate_charges_with_resp(self, mol_name: str, pdb_file: str, resp_path: str) -> None:
        """使用RESP2计算分子的电荷"""
        
        # 硬编码检查常见已知离子
        known_ions = ["Li", "Na", "K", "PF6", "TFSI", "BF4", "FSI", "Li+", "Na+", "K+", "PF6-", "TFSI-", "BF4-", "FSI-"]
        if mol_name in known_ions:
            self.logger.info(f"跳过已知离子 {mol_name} 的电荷生成")
            return
            
        # 检查配置中的离子清单
        for cat in self.config.get("cations", []):
            if cat.get("name") == mol_name:
                self.logger.info(f"跳过阳离子 {mol_name} 的电荷生成")
                return
        
        for an in self.config.get("anions", []):
            if an.get("name") == mol_name:
                self.logger.info(f"跳过阴离子 {mol_name} 的电荷生成") 
                return
        
        self.logger.info(f"为溶剂分子 {mol_name} 使用RESP2生成电荷")
        
        try:
            # 确保PDB文件存在
            if not os.path.exists(pdb_file):
                self.logger.error(f"PDB文件不存在: {pdb_file}")
                return
                
            # 确保RESP目录存在
            if not os.path.exists(resp_path):
                self.logger.error(f"RESP目录不存在: {resp_path}")
                return
            
            # 确保电荷保存目录存在
            charge_save_dir = self.charge_save_path
            if not charge_save_dir:
                charge_save_dir = os.path.join(os.path.dirname(pdb_file), "charges")
            os.makedirs(charge_save_dir, exist_ok=True)
            
            # 复制PDB文件到RESP目录
            dest_pdb = os.path.join(resp_path, f"{mol_name}.pdb")
            shutil.copy2(pdb_file, dest_pdb)
            self.logger.info(f"已复制PDB文件到RESP目录: {pdb_file} -> {dest_pdb}")
            
            # 运行RESP脚本
            resp_script = os.path.join(resp_path, "RESP2.sh")
            if os.path.exists(resp_script):
                # 切换到RESP目录
                original_dir = os.getcwd()
                os.chdir(resp_path)
                
                # 先修改RESP2.sh文件，直接替换Gaussian变量为绝对路径
                g16_path = "/public/software/g16/"
                
                # 创建RESP2.sh的备份
                resp_script_backup = os.path.join(resp_path, "RESP2.sh.bak")
                shutil.copy2(resp_script, resp_script_backup)
                
                # 读取RESP2.sh内容
                with open(resp_script, 'r') as f:
                    script_content = f.read()
                
                # 替换Gaussian=g16行为绝对路径
                script_content = script_content.replace('Gaussian=g16', f'Gaussian="{g16_path}"')
                script_content = script_content.replace('Gaussian=g09', f'Gaussian="{g16_path}"')
                
                # 写回修改后的内容
                with open(resp_script, 'w') as f:
                    f.write(script_content)
                
                self.logger.info(f"已修改RESP2.sh，将Gaussian变量设置为绝对路径: {g16_path}")
                
                # 使用绝对路径运行命令
                cmd = f"""
export g16root=/public/software/g16
export GAUSS_EXEDIR=/public/software/g16
export PATH=$GAUSS_EXEDIR:$PATH
export Multiwfnpath=/public/home/xiaoji/software/Multiwfn_3.7_bin_Linux
export PATH=$Multiwfnpath:$PATH
bash RESP2.sh {mol_name}.pdb
"""
                self.logger.info(f"执行RESP命令: {cmd}")
                
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                stdout, stderr = process.communicate()
                
                # 恢复原始目录
                os.chdir(original_dir)
                
                # 检查命令执行是否成功
                if process.returncode != 0:
                    self.logger.error(f"RESP命令执行失败: {stderr}")
                    return  # 直接返回，不生成进一步的警告
                
                self.logger.info(f"RESP命令执行成功: {stdout}")
                
                # 恢复原始RESP2.sh
                shutil.copy2(resp_script_backup, resp_script)
                
                # 检查是否生成了电荷文件 (.chg文件)
                charge_file = os.path.join(resp_path, f"{mol_name}.chg")
                if os.path.exists(charge_file):
                    # 复制电荷文件到电荷保存目录
                    dest_charge = os.path.join(charge_save_dir, f"{mol_name}.chg")
                    shutil.copy2(charge_file, dest_charge)
                    self.logger.info(f"已保存电荷文件: {dest_charge}")
                    
                    # 复制PDB文件到电荷保存目录
                    dest_pdb_save = os.path.join(charge_save_dir, f"{mol_name}.pdb")
                    shutil.copy2(pdb_file, dest_pdb_save)
                    self.logger.info(f"已保存PDB文件: {dest_pdb_save}")
                    
                    # 立即更新exist_molecules.json
                    self._update_exist_solvent_json(mol_name)
                else:
                    self.logger.warning(f"RESP未生成电荷文件: {charge_file}")
            else:
                self.logger.error(f"找不到RESP脚本: {resp_script}")
        
        except Exception as e:
            self.logger.error(f"生成RESP电荷时出错: {str(e)}")
            # 记录更详细的异常信息
            import traceback
            self.logger.error(f"详细错误信息: {traceback.format_exc()}")

    def _update_exist_solvent_json(self, mol_name: str) -> None:
        """更新exist_solvent.json文件，添加新的溶剂信息
        
        Args:
            mol_name: 分子名称
        """
        try:
            # 获取分子的SMILE字符串
            smile = ""
            for solvent in self.config.get("solvents", []):
                if solvent.get("name") == mol_name:
                    smile = solvent.get("smile", "")
                    break
            
            if not smile:
                self.logger.warning(f"未找到分子 {mol_name} 的SMILE字符串")
                return
            
            # 读取现有的exist_solvent.json
            existing_molecules_file = os.path.join(self.resp_path, "existing_molecules.json")
            self.logger.info(f"使用existing_molecules.json文件路径: {existing_molecules_file}")
                
            # 确保目录存在
            os.makedirs(os.path.dirname(existing_molecules_file), exist_ok=True)
            
            # 读取现有数据
            existing_data = {
                "solvent": {},
                "cation": {},
                "anion": {}
            }
            if os.path.exists(existing_molecules_file):
                with open(existing_molecules_file, 'r') as f:
                    try:
                        existing_data = json.load(f)
                    except json.JSONDecodeError:
                        self.logger.warning("existing_molecules.json格式错误，将创建新文件")
                        existing_data = {
                            "solvent": {},
                            "cation": {},
                            "anion": {}
                        }
            
            # 设置电荷文件路径(.chg而不是.charge)
            charge_file_path = os.path.join(self.charge_save_path, f"{mol_name}.chg")
            if not os.path.exists(charge_file_path):
                # 检查RESP目录
                resp_charge_path = os.path.join(self.resp_path, f"{mol_name}.chg")
                if os.path.exists(resp_charge_path):
                    charge_file_path = resp_charge_path
                else:
                    self.logger.warning(f"找不到 {mol_name} 的.chg电荷文件")
                    charge_file_path = f"{mol_name}.chg"  # 使用相对路径
            
            # 更新数据 - 添加到溶剂部分
            existing_data["solvent"][mol_name] = {
                "name": mol_name,
                "smile": smile,
                "charge_file": charge_file_path,
                "last_updated": datetime.now().isoformat()
            }
            
            # 保存更新后的数据
            with open(existing_molecules_file, 'w') as f:
                json.dump(existing_data, f, indent=2)
            
            self.logger.info(f"已更新existing_molecules.json，添加溶剂分子: {mol_name}，电荷文件: {charge_file_path}")
            
        except Exception as e:
            self.logger.error(f"更新existing_molecules.json时出错: {e}")
            import traceback
            self.logger.error(f"详细错误信息: {traceback.format_exc()}")

    def create_lammps_files(self):
        """创建LAMMPS输入文件"""
        try:
            # 确保输出目录存在
            if not self.output_dir:
                self.output_dir = os.path.join(self.work_dir, "output") if self.work_dir else "output"
            os.makedirs(self.output_dir, exist_ok=True)
            
            # 记录当前目录，以便后续恢复
            current_dir = os.getcwd()
            
            # 切换到输出目录
            os.chdir(self.output_dir)
            self.logger.info(f"切换到输出目录: {self.output_dir}")
            
            # 调用已有方法或创建临时方法来生成分子文件
            if hasattr(self, '_create_molecule_files') and callable(getattr(self, '_create_molecule_files')):
                self.logger.info("调用_create_molecule_files方法")
                self._create_molecule_files()
            else:
                self.logger.warning("_create_molecule_files方法未定义，尝试使用generate_input_files")
                if hasattr(self, 'generate_input_files') and callable(getattr(self, 'generate_input_files')):
                    self.logger.info("调用generate_input_files方法")
                    result = self.generate_input_files()
                    if result:
                        self.logger.info(f"输入文件生成成功: {result}")
                        
                        # 更新生成的文件列表
                        if not hasattr(self, 'generated_files') or self.generated_files is None:
                            self.generated_files = {}
                        if isinstance(result, dict):
                            self.generated_files.update(result)
                        
                        # 已经生成了所有文件，可以直接返回
                        self.lammps_files_created = True
                        
                        # 恢复原始目录
                        os.chdir(current_dir)
                        self.logger.info(f"恢复原始目录: {current_dir}")
                        
                        # 更新关键文件路径
                        if 'system_pdb' in result and result['system_pdb'] and os.path.exists(result['system_pdb']):
                            self.pdb_file = result['system_pdb']
                        if 'lammps_data' in result and result['lammps_data'] and os.path.exists(result['lammps_data']):
                            self.data_file = result['lammps_data']
                        if 'lammps_input' in result and result['lammps_input'] and os.path.exists(result['lammps_input']):
                            self.in_file = result['lammps_input']
                        
                        return self.generated_files
                else:
                    self.logger.error("无法找到方法来生成分子文件")
                    raise ValueError("无法找到方法来生成分子文件")
            
            # 创建系统LT文件
            system_lt_path = os.path.join(self.output_dir, "system.lt")
            if hasattr(self, '_create_system_lt_file') and callable(getattr(self, '_create_system_lt_file')):
                self.logger.info("调用_create_system_lt_file方法")
                self._create_system_lt_file(system_lt_path)
            else:
                self.logger.warning("_create_system_lt_file方法未定义，将跳过")
            
            # 使用moltemplate生成LAMMPS数据文件
            if hasattr(self, '_run_moltemplate') and callable(getattr(self, '_run_moltemplate')):
                self.logger.info("调用_run_moltemplate方法")
                self._run_moltemplate()
            else:
                self.logger.warning("_run_moltemplate方法未定义，将跳过")
            
            # 检查生成的文件
            data_file = os.path.join(self.output_dir, "system.data")
            in_file = os.path.join(self.output_dir, "system.in")
            in_settings_file = os.path.join(self.output_dir, "system.in.settings")
            in_init_file = os.path.join(self.output_dir, "system.in.init")
            pdb_file = os.path.join(self.output_dir, "system.pdb")
            
            # 记录生成的文件路径
            self.generated_files = {
                "lammps_data": data_file if os.path.exists(data_file) else None,
                "lammps_input": in_file if os.path.exists(in_file) else None,
                "system_settings": in_settings_file if os.path.exists(in_settings_file) else None,
                "system_init": in_init_file if os.path.exists(in_init_file) else None,
                "system_pdb": pdb_file if os.path.exists(pdb_file) else None,
                "system_lt": system_lt_path if os.path.exists(system_lt_path) else None
            }
            
            # 记录关键文件路径
            self.data_file = data_file if os.path.exists(data_file) else None
            self.in_file = in_file if os.path.exists(in_file) else None
            self.pdb_file = pdb_file if os.path.exists(pdb_file) else None
            
            # 记录文件生成状态
            for name, path in self.generated_files.items():
                if path and os.path.exists(path):
                    self.logger.info(f"生成的{name}文件: {path}")
                else:
                    self.logger.warning(f"{name}文件不存在或路径为None")
            
            # 创建LAMMPS运行脚本
            if hasattr(self, '_create_run_script') and callable(getattr(self, '_create_run_script')):
                self.logger.info("调用_create_run_script方法")
                self._create_run_script()
            else:
                self.logger.warning("_create_run_script方法未定义，将跳过")
            
            # 创建作业提交脚本
            if hasattr(self, '_create_job_script') and callable(getattr(self, '_create_job_script')):
                self.logger.info("调用_create_job_script方法")
                self._create_job_script()
            else:
                self.logger.warning("_create_job_script方法未定义，将跳过")
            
            # 恢复原始目录
            os.chdir(current_dir)
            self.logger.info(f"恢复原始目录: {current_dir}")
            
            # 标记LAMMPS文件已创建
            self.lammps_files_created = True
            
            self.logger.info("输入文件生成完成")
            return self.generated_files
        
        except Exception as e:
            self.logger.error(f"创建LAMMPS输入文件失败: {str(e)}")
            self.logger.error(traceback.format_exc())
            # 恢复原始目录
            if 'current_dir' in locals():
                os.chdir(current_dir)
                self.logger.info(f"恢复原始目录: {current_dir}")
            raise ValueError(f"创建LAMMPS输入文件失败: {str(e)}")

    def submit_job(self, callback=None):
        """提交SLURM作业
        
        Args:
            callback: 可选的回调函数，用于在获取作业ID后执行
            
        Returns:
            str: SLURM作业ID
        """
        try:
            self.logger.info("开始提交SLURM作业")
            
            # 确保作业脚本存在
            job_script = os.path.join(self.output_dir, "job.sh")
            if not os.path.exists(job_script):
                self.logger.warning(f"作业脚本不存在: {job_script}")
                
                # 尝试创建作业脚本
                if hasattr(self, '_create_job_script') and callable(getattr(self, '_create_job_script')):
                    self.logger.info("尝试创建作业脚本")
                    self._create_job_script()
                    
                    # 再次检查
                    if not os.path.exists(job_script):
                        error_msg = f"无法创建作业脚本: {job_script}"
                        self.logger.error(error_msg)
                        raise ValueError(error_msg)
                else:
                    error_msg = f"作业脚本不存在且无法创建: {job_script}"
                    self.logger.error(error_msg)
                    raise ValueError(error_msg)
            
            # 保存当前目录
            current_dir = os.getcwd()
            
            # 切换到包含作业脚本的目录
            script_dir = os.path.dirname(job_script)
            os.chdir(script_dir)
            self.logger.info(f"切换到作业脚本目录: {script_dir}")
            
            # 获取作业脚本文件名（不包含路径）
            script_name = os.path.basename(job_script)
            
            # 使用sbatch提交作业
            self.logger.info(f"执行命令: sbatch {script_name}")
            try:
                result = subprocess.run(["sbatch", script_name], capture_output=True, text=True)
            except Exception as e:
                self.logger.error(f"执行sbatch命令失败: {str(e)}")
                # 恢复原始目录
                os.chdir(current_dir)
                raise ValueError(f"执行sbatch命令失败: {str(e)}")
            
            # 恢复原始目录
            os.chdir(current_dir)
            self.logger.info(f"恢复原始目录: {current_dir}")
            
            # 检查提交结果
            if result.returncode == 0:
                # 从输出中提取作业ID
                output = result.stdout
                self.logger.info(f"sbatch输出: {output}")
                
                # 使用正则表达式提取作业ID
                match = re.search(r"Submitted batch job (\d+)", output)
                if match:
                    job_id = match.group(1)
                    self.logger.info(f"成功提交SLURM作业: {job_id}")
                    
                    # 保存作业ID
                    self.slurm_job_id = job_id
                    
                    # 如果提供了回调函数，执行它
                    if callback and callable(callback):
                        try:
                            self.logger.info(f"执行回调函数")
                            callback(job_id)
                            self.logger.info(f"回调函数执行成功")
                        except Exception as e:
                            self.logger.error(f"执行回调函数时出错: {str(e)}")
                            # 继续执行，不要因回调失败而终止整个过程
                    
                    return job_id
                else:
                    error_msg = f"无法从sbatch输出中提取作业ID: {output}"
                    self.logger.error(error_msg)
                    raise ValueError(error_msg)
            else:
                error_msg = f"提交作业失败: {result.stderr}"
                self.logger.error(error_msg)
                raise ValueError(error_msg)
        
        except Exception as e:
            self.logger.error(f"提交SLURM作业时出错: {str(e)}")
            self.logger.error(traceback.format_exc())
            # 返回None而不是抛出异常，这样即使作业提交失败，工作流也能继续执行
            return None

# 提供便捷的函数用于从INP文件运行工作流
def run_from_inp_file(inp_file_path: str, work_dir: str = None, output_dir: str = None, 
                   user_id: str = None, project_id: str = None, formulation_id: str = None) -> Dict[str, Any]:
    """从INP文件运行完整的电解液计算工作流
    
    Args:
        inp_file_path: 输入文件路径
        work_dir: 工作目录
        output_dir: 输出目录
        user_id: 用户ID
        project_id: 项目ID
        formulation_id: 配方ID
        
    Returns:
        Dict: 计算结果
    """
    try:
        # 创建工作流实例
        workflow = MolyteWorkflow(
            inp_file_path=inp_file_path,
            work_dir=work_dir,
            output_dir=output_dir,
            user_id=user_id,
            project_id=project_id,
            formulation_id=formulation_id
        )
        
        # 运行完整工作流
        result = workflow.run()
        return result
        
    except Exception as e:
        # 记录错误并返回
        logging.error(f"运行工作流时出错: {str(e)}", exc_info=True)
        return {
            'status': 'failed',
            'error': str(e)
        }


# 提供便捷的函数用于从INP内容字符串运行工作流
def run_from_inp_content(inp_content: str, 
                        work_dir: Optional[str] = None,
                        output_dir: Optional[str] = None,
                        user_id: Optional[str] = None,
                        project_id: Optional[str] = None,
                        formulation_id: Optional[str] = None) -> Dict[str, Any]:
    """从INP内容字符串运行电解液工作流
    
    Args:
        inp_content: INP内容字符串
        work_dir: 工作目录
        output_dir: 输出目录
        user_id: 用户ID
        project_id: 项目ID
        formulation_id: 配方ID
        
    Returns:
        Dict: 计算结果
    """
    workflow = MolyteWorkflow(
        inp_content=inp_content,
        work_dir=work_dir,
        output_dir=output_dir,
        user_id=user_id,
        project_id=project_id,
        formulation_id=formulation_id
    )
    return workflow.run() 