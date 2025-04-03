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
        
        # 确保必要的目录存在
        os.makedirs(self.initial_salts_path, exist_ok=True)
        os.makedirs(self.initial_solvent_path, exist_ok=True)
        
        # 是否使用RESP电荷
        self.use_resp_charges = paths.get('use_resp_charges', True)
        
        # 记录路径信息
        self.logger.info("外部程序路径设置:")
        self.logger.info(f"  ligpargen: {self.ligpargen_path}")
        self.logger.info(f"  packmol: {self.packmol_path}")
        self.logger.info(f"  ltemplify: {self.ltemplify_path}")
        self.logger.info(f"  moltemplate: {self.moltemplate_path}")
        self.logger.info(f"  resp: {self.resp_path}")
        self.logger.info(f"  通用数据路径: {self.common_data_path}")
        self.logger.info(f"  电荷保存路径: {self.charge_save_path}")
        self.logger.info(f"  初始盐路径: {self.initial_salts_path}")
        self.logger.info(f"  初始溶剂路径: {self.initial_solvent_path}")
    
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
        
        # 创建分子文件目录
        molecule_dir = os.path.join(self.work_dir, "molecules")
        os.makedirs(molecule_dir, exist_ok=True)
        
        # 处理所有分子
        self.logger.info("处理分子文件")
        molecules_info = []
        
        # 合并相同阳离子
        self.logger.info("合并相同名称的阳离子")
        cations_map = {}
        for cation in self.formulation_data.get("cations", []):
            molecule_name = cation.get("name")
            if not molecule_name:
                self.logger.warning("发现没有名称的阳离子，跳过")
                continue
                
            if molecule_name in cations_map:
                existing_cation = cations_map[molecule_name]
                existing_number = existing_cation.get("number", 0)
                new_number = cation.get("number", 0)
                existing_cation["number"] = existing_number + new_number
                self.logger.info(f"合并阳离子 {molecule_name}: {existing_number} + {new_number} = {existing_cation['number']}")
            else:
                cations_map[molecule_name] = cation.copy()
                self.logger.info(f"添加新阳离子 {molecule_name}: {cation.get('number', 0)}")
        
        # 合并相同阴离子
        self.logger.info("合并相同名称的阴离子")
        anions_map = {}
        for anion in self.formulation_data.get("anions", []):
            molecule_name = anion.get("name")
            if not molecule_name:
                self.logger.warning("发现没有名称的阴离子，跳过")
                continue
                
            if molecule_name in anions_map:
                existing_anion = anions_map[molecule_name]
                existing_number = existing_anion.get("number", 0)
                new_number = anion.get("number", 0)
                existing_anion["number"] = existing_number + new_number
                self.logger.info(f"合并阴离子 {molecule_name}: {existing_number} + {new_number} = {existing_anion['number']}")
            else:
                anions_map[molecule_name] = anion.copy()
                self.logger.info(f"添加新阴离子 {molecule_name}: {anion.get('number', 0)}")
        
        # 处理合并后的阳离子
        merged_cations = list(cations_map.values())
        self.logger.info(f"处理合并后的阳离子: {len(merged_cations)}个")
        for cation in merged_cations:
            molecule_name = cation.get("name")
            self.logger.info(f"处理阳离子: {molecule_name}, 数量: {cation.get('number', 0)}")
            # 添加类型标识
            cation["type"] = "cation"
            try:
                files = self._generate_molecule_files(cation, molecule_dir)
                if files:
                    sanitized_name = molecule_name.replace('+', '').replace('-', '')
                    cation["sanitized_name"] = sanitized_name
                    cation["files"] = files
                    molecules_info.append(cation)
                    self.logger.info(f"阳离子 {molecule_name} 文件处理完成")
            except Exception as e:
                self.logger.error(f"阳离子 {molecule_name} 文件生成失败: {str(e)}")
                raise RuntimeError(f"阳离子 {molecule_name} 文件生成失败: {str(e)}")
                
        # 处理合并后的阴离子
        merged_anions = list(anions_map.values())
        self.logger.info(f"处理合并后的阴离子: {len(merged_anions)}个")
        for anion in merged_anions:
            molecule_name = anion.get("name")
            self.logger.info(f"处理阴离子: {molecule_name}, 数量: {anion.get('number', 0)}")
            # 添加类型标识
            anion["type"] = "anion"
            try:
                files = self._generate_molecule_files(anion, molecule_dir)
                if files:
                    sanitized_name = molecule_name.replace('+', '').replace('-', '')
                    anion["sanitized_name"] = sanitized_name
                    anion["files"] = files
                    molecules_info.append(anion)
                    self.logger.info(f"阴离子 {molecule_name} 文件处理完成")
            except Exception as e:
                self.logger.error(f"阴离子 {molecule_name} 文件生成失败: {str(e)}")
                raise RuntimeError(f"阴离子 {molecule_name} 文件生成失败: {str(e)}")
                
        # 处理溶剂
        for solvent in self.formulation_data.get("solvents", []):
            molecule_name = solvent.get("name")
            self.logger.info(f"处理溶剂: {molecule_name}")
            # 添加类型标识
            solvent["type"] = "solvent"
            try:
                files = self._generate_molecule_files(solvent, molecule_dir)
                if files:
                    sanitized_name = molecule_name.replace('+', '').replace('-', '')
                    solvent["sanitized_name"] = sanitized_name
                    solvent["files"] = files
                    molecules_info.append(solvent)
                    self.logger.info(f"溶剂 {molecule_name} 文件处理完成")
            except Exception as e:
                self.logger.error(f"溶剂 {molecule_name} 文件生成失败: {str(e)}")
                raise RuntimeError(f"溶剂 {molecule_name} 文件生成失败: {str(e)}")
        
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
    
    def _generate_molecule_files(self, molecule, molecule_dir):
        """处理分子文件，对于阴阳离子直接检查initial_salts中，对于溶剂检查initial_solvent或用ligpargen生成"""
        
        # 获取分子名称和信息
        molecule_name = molecule.get("name")
        smile = molecule.get("smile", "")
        
        # 根据调用的情况确定分子类型，不再尝试从名称判断
        # 检查molecule对象来自哪个列表
        molecule_type = molecule.get("type", "unknown")
        
        # 为阴阳离子检查initial_salts目录
        if molecule_type == "cation" or molecule_type == "anion":
            self.logger.info(f"为{molecule_type} {molecule_name} 检查initial_salts目录")
            salt_files = {}
            
            # 检查必要的文件扩展名: lt 和 pdb
            for ext in ['lt', 'pdb']:
                file_path = os.path.join(self.initial_salts_path, f"{molecule_name}.{ext}")
                if os.path.exists(file_path):
                    # 复制到目标目录
                    dst_file = os.path.join(molecule_dir, f"{molecule_name}.{ext}")
                    shutil.copy2(file_path, dst_file)
                    salt_files[ext] = dst_file
                    self.logger.info(f"已复制文件: {file_path} -> {dst_file}")
                else:
                    self.logger.error(f"找不到{molecule_name}的{ext}文件: {file_path}")
            
            # 验证是否找到了所有必要的文件
            if 'lt' in salt_files and 'pdb' in salt_files:
                self.logger.info(f"成功找到{molecule_type} {molecule_name} 的所有必要文件")
                return salt_files
            else:
                missing_files = []
                if 'lt' not in salt_files:
                    missing_files.append('lt')
                if 'pdb' not in salt_files:
                    missing_files.append('pdb')
                
                self.logger.error(f"{molecule_type} {molecule_name} 缺少必要的文件: {', '.join(missing_files)}")
                raise RuntimeError(f"{molecule_type} {molecule_name} 缺少必要的文件: {', '.join(missing_files)}")
        
        # 对于溶剂，先检查initial_solvent目录
        elif molecule_type == "solvent":
            self.logger.info(f"处理溶剂: {molecule_name}")
            initial_solvent_dir = self.initial_solvent_path
            
            if hasattr(self, 'initial_solvent_path') and os.path.exists(initial_solvent_dir):
                # 检查initial_solvent目录中是否有该溶剂
                self.logger.info(f"检查initial_solvent目录中是否有溶剂 {molecule_name}")
                solvent_files = {}
                
                for ext in ['pdb', 'lt', 'lmp', 'chg']:
                    file_path = os.path.join(initial_solvent_dir, f"{molecule_name}.{ext}")
                    if os.path.exists(file_path):
                        # 复制到目标目录
                        dst_file = os.path.join(molecule_dir, f"{molecule_name}.{ext}")
                        shutil.copy2(file_path, dst_file)
                        solvent_files[ext] = dst_file
                        self.logger.info(f"从initial_solvent复制文件: {file_path} -> {dst_file}")
                
                # 如果找到了必要的文件，直接返回
                if 'pdb' in solvent_files and 'lt' in solvent_files:
                    self.logger.info(f"在initial_solvent中成功找到溶剂 {molecule_name} 的必要文件")
                    # 保存smile以便后续查找
                    if smile:
                        self._update_exist_solvent_json(molecule_name)
                    return solvent_files
            
            # 如果没找到溶剂或缺少必要文件，使用ligpargen生成
            if not smile:
                self.logger.error(f"无法生成溶剂 {molecule_name}：未提供SMILE字符串")
                raise RuntimeError(f"无法生成溶剂 {molecule_name}：未提供SMILE字符串")
            
            # 创建临时目录用于生成分子文件
            temp_dir = os.path.join(molecule_dir, "temp_gen")
            os.makedirs(temp_dir, exist_ok=True)
            
            # 使用ligpargen处理SMILE
            self.logger.info(f"使用ligpargen处理SMILE生成溶剂: {smile}")
            try:
                # 生成溶剂文件
                result = self._generate_molecule_files_with_ligpargen(
                    {'name': molecule_name, 'smile': smile}, 
                    temp_dir
                )
                
                if not result:
                    raise RuntimeError(f"使用ligpargen生成溶剂 {molecule_name} 失败")
                
                # 移动生成的文件到目标目录
                for ext, file_path in result.items():
                    dst_file = os.path.join(molecule_dir, f"{molecule_name}.{ext}")
                    shutil.copy2(file_path, dst_file)
                    result[ext] = dst_file
                    self.logger.info(f"复制生成的文件: {file_path} -> {dst_file}")
                
                # 如果有lmp文件但没有lt文件，从lmp生成lt
                if 'lmp' in result and 'lt' not in result:
                    lt_file = os.path.join(molecule_dir, f"{molecule_name}.lt")
                    self._generate_lt_from_lmp(result['lmp'], lt_file, molecule_name)
                    if os.path.exists(lt_file):
                        result['lt'] = lt_file
                        self.logger.info(f"从LMP生成LT文件: {lt_file}")
                
                # 为resp生成的溶剂保存到initial_solvent
                if 'pdb' in result and 'lt' in result:
                    self.logger.info(f"将生成的溶剂 {molecule_name} 保存到initial_solvent")
                    
                    # 确保initial_solvent目录存在
                    os.makedirs(initial_solvent_dir, exist_ok=True)
                    
                    # 复制文件到initial_solvent
                    for ext, file_path in result.items():
                        dst_file = os.path.join(initial_solvent_dir, f"{molecule_name}.{ext}")
                        shutil.copy2(file_path, dst_file)
                        self.logger.info(f"复制到initial_solvent: {file_path} -> {dst_file}")
                    
                    # 保存smile信息
                    if smile:
                        self._update_exist_solvent_json(molecule_name)
                
                # 清理临时目录
                shutil.rmtree(temp_dir, ignore_errors=True)
                
                return result
                
            except Exception as e:
                self.logger.error(f"生成溶剂 {molecule_name} 时出错: {str(e)}")
                raise RuntimeError(f"生成溶剂 {molecule_name} 时出错: {str(e)}")
        
        # 未知分子类型
        else:
            self.logger.error(f"未知分子类型: {molecule_type} - {molecule_name}")
            raise RuntimeError(f"未知分子类型: {molecule_type} - {molecule_name}")

    def _get_molecule_files(self, directory):
        """扫描目录并返回找到的分子文件
        
        Args:
            directory: 要扫描的目录路径
            
        Returns:
            Dict: 文件类型到文件路径的映射
        """
        self.logger.info(f"扫描目录中的分子文件: {directory}")
        
        if not os.path.exists(directory):
            self.logger.error(f"目录不存在: {directory}")
            return {}
            
        result_files = {}
        
        # 常见的分子文件扩展名
        extensions = ["pdb", "lt", "chg", "lmp"]
        
        # 扫描目录中的文件
        for file in os.listdir(directory):
            file_path = os.path.join(directory, file)
            if os.path.isfile(file_path):
                # 获取文件扩展名
                _, ext = os.path.splitext(file)
                ext = ext.lower().lstrip('.')
                
                if ext in extensions:
                    result_files[ext] = file_path
                    self.logger.info(f"找到分子文件: {ext} => {file_path}")
        
        return result_files

    def generate_molecule_from_smile(self, smile, molecule_name, output_dir):
        """从SMILE字符串生成分子文件
        
        Args:
            smile: 分子的SMILE字符串
            molecule_name: 分子名称
            output_dir: 输出目录路径
            
        Returns:
            bool: 成功返回True，失败返回False
        """
        self.logger.info(f"从SMILE生成分子文件: {molecule_name}, SMILE: {smile}")
        
        # 创建分子信息字典
        molecule_info = {
            'name': molecule_name,
            'smile': smile
        }
        
        # 调用LigParGen处理函数
        result_files = self._generate_molecule_files_with_ligpargen(molecule_info, output_dir)
        
        # 检查是否成功生成文件
        if result_files and len(result_files) > 0:
            self.logger.info(f"成功从SMILE生成分子文件: {molecule_name}, 文件: {result_files}")
            return True
        else:
            self.logger.error(f"从SMILE生成分子文件失败: {molecule_name}")
            return False

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
        """使用LigParGen生成溶剂分子文件，并使用resp2生成电荷
        
        Args:
            molecule_info: 分子信息字典，包含name和smile
            molecule_dir: 输出分子目录
            
        Returns:
            Dict: 包含生成文件路径的字典，如果失败则返回空字典
        """
        # 打印版本标识，帮助验证使用了哪个版本的代码
        self.logger.info("=== 使用修改后的LigParGen处理函数 ===")
        
        molecule_name = molecule_info.get('name', '')
        molecule_smile = molecule_info.get('smile', '')
        
        if not molecule_name or not molecule_smile:
            self.logger.warning(f"分子名称或SMILE字符串为空，无法使用LigParGen生成: {molecule_name}")
            return {}
        
        # 创建分子子目录
        sanitized_name = molecule_name.replace('+', '').replace('-', '')
        molecule_subdir = os.path.join(molecule_dir, sanitized_name)
        os.makedirs(molecule_subdir, exist_ok=True)
        
        try:
            # 切换到输出目录
            original_dir = os.getcwd()
            os.chdir(molecule_dir)
            
            # 准备输出文件路径
            pdb_file = os.path.join(molecule_subdir, f"{sanitized_name}.pdb")
            lmp_file = os.path.join(molecule_subdir, f"{sanitized_name}.lmp")
            lt_file = os.path.join(molecule_subdir, f"{sanitized_name}.lt")
            
            # 构建LigParGen命令
            cmd_str = f"{self.ligpargen_path} -s '{molecule_smile}' -n {sanitized_name} -r MOL -c 0 -o 0 -cgen CM1A"
            
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
            self.logger.info(f"当前目录文件列表: {current_files}")
            
            # 查找PDB文件
            pdb_files = [f for f in current_files if f.endswith('.pdb')]
            if pdb_files:
                # 使用找到的第一个.pdb文件
                src_pdb = os.path.join(os.getcwd(), pdb_files[0])
                self.logger.info(f"使用PDB文件: {src_pdb}")
                shutil.copy2(src_pdb, pdb_file)
                result_files['pdb'] = pdb_file
            else:
                self.logger.error(f"未找到任何PDB文件")
                return {}
            
            # 查找LMP文件
            lmp_files = [f for f in current_files if f.endswith('.lmp')]
            if lmp_files:
                # 使用找到的第一个.lmp文件
                src_lmp = os.path.join(os.getcwd(), lmp_files[0])
                self.logger.info(f"使用LMP文件: {src_lmp}")
                shutil.copy2(src_lmp, lmp_file)
                result_files['lmp'] = lmp_file
            else:
                self.logger.error(f"未找到任何LMP文件")
                return {}
            
            # 使用RESP2生成电荷
            charge_file = os.path.join(molecule_subdir, f"{sanitized_name}.chg")
            if os.path.exists(self.resp_path):
                self.logger.info(f"使用RESP2为分子{sanitized_name}生成电荷文件")
                try:
                    # 调用RESP生成电荷
                    self._generate_charges_with_resp(sanitized_name, result_files['pdb'], self.resp_path)
                    
                    # 检查是否生成了电荷文件
                    resp_charge_file = os.path.join(self.resp_path, f"{sanitized_name}.chg")
                    if os.path.exists(resp_charge_file):
                        self.logger.info(f"RESP生成了电荷文件: {resp_charge_file}")
                        # 复制并重命名为.chg文件
                        shutil.copy2(resp_charge_file, charge_file)
                        result_files['chg'] = charge_file
                        
                        # 更新LMP文件中的电荷
                        if 'lmp' in result_files:
                            try:
                                self.logger.info(f"使用RESP电荷更新LMP文件中的电荷")
                                replace_charges_in_lmp(result_files['lmp'], charge_file)
                                self.logger.info(f"成功更新LMP文件中的电荷")
                            except Exception as e:
                                self.logger.error(f"更新LMP文件电荷时出错: {str(e)}")
                except Exception as e:
                    self.logger.error(f"使用RESP生成电荷时出错: {str(e)}")
            
            # 从LMP生成LT文件
            if 'lmp' in result_files:
                try:
                    self.logger.info(f"从LMP文件生成LT文件")
                    lt_file_path = self._generate_lt_from_lmp(result_files['lmp'], lt_file, sanitized_name)
                    if os.path.exists(lt_file_path):
                        result_files['lt'] = lt_file_path
                        self.logger.info(f"成功从LMP生成LT文件: {lt_file_path}")
                except Exception as e:
                    self.logger.error(f"从LMP生成LT文件失败: {e}")
            
            # 返回到原始目录
            os.chdir(original_dir)
            
            # 检查必要的文件是否生成
            if 'pdb' not in result_files or 'lmp' not in result_files:
                self.logger.error(f"缺少必要的文件: PDB或LMP")
                return {}
            
            return result_files
            
        except Exception as e:
            # 返回到原始目录
            if original_dir:
                os.chdir(original_dir)
            self.logger.error(f"使用LigParGen生成分子文件失败: {str(e)}")
            return {}

    def _generate_lt_from_lmp(self, lmp_file, lt_file, molecule_name):
        """使用ltemplify从LMP文件生成LT文件，并确保LT文件中的电荷与CHG文件一致
        
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
                
                # 寻找相应的CHG文件，确保LT文件中的电荷与CHG文件一致
                chg_file = None
                
                # 先检查同一目录下是否有对应的chg文件
                possible_chg_file = os.path.join(os.path.dirname(lmp_file), f"{molecule_name}.chg")
                if os.path.exists(possible_chg_file):
                    chg_file = possible_chg_file
                    self.logger.info(f"找到相应的电荷文件: {chg_file}")
                
                # 如果没找到，检查LMP文件目录的父目录
                if not chg_file:
                    parent_dir = os.path.dirname(os.path.dirname(lmp_file))
                    possible_chg_file = os.path.join(parent_dir, f"{molecule_name}.chg")
                    if os.path.exists(possible_chg_file):
                        chg_file = possible_chg_file
                        self.logger.info(f"在父目录找到相应的电荷文件: {chg_file}")
                
                # 如果找到了电荷文件，则更新LT文件中的电荷
                if chg_file:
                    self.logger.info(f"更新LT文件中的电荷，使用电荷文件: {chg_file}")
                    self._update_charges_in_lt(lt_file, chg_file)
                else:
                    self.logger.warning(f"未找到相应的电荷文件，LT文件中的电荷可能不准确")
                
                self.logger.info(f"成功生成LT文件: {lt_file}")
                return lt_file
            else:
                error_msg = f"ltemplify未能生成LT文件: {lt_file}"
                self.logger.error(error_msg)
                raise FileNotFoundError(error_msg)
                
        except Exception as e:
            self.logger.error(f"执行ltemplify时出错: {e}")
            raise

    def _update_charges_in_lt(self, lt_file, chg_file):
        """更新LT文件中的电荷信息，使其与CHG文件一致
        
        Args:
            lt_file: LT文件路径
            chg_file: 电荷文件路径
        """
        self.logger.info(f"更新LT文件中的电荷: {lt_file}, 使用电荷文件: {chg_file}")
        
        try:
            # 读取CHG文件中的电荷
            charges = []
            with open(chg_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split()
                        if len(parts) > 0:
                            try:
                                charge = float(parts[-1])
                                charges.append(charge)
                            except ValueError:
                                pass
            
            if not charges:
                self.logger.warning(f"从电荷文件未能读取到有效的电荷值: {chg_file}")
                return
                
            self.logger.info(f"从电荷文件读取到 {len(charges)} 个电荷值")
            
            # 读取LT文件内容
            with open(lt_file, 'r') as f:
                lt_content = f.readlines()
            
            # 查找write("Data Atoms")部分，这是定义原子电荷的地方
            atom_section_start = None
            atom_section_end = None
            
            for i, line in enumerate(lt_content):
                if 'write("Data Atoms")' in line:
                    atom_section_start = i + 1
                    break
            
            if atom_section_start is None:
                self.logger.warning(f"在LT文件中未找到原子定义部分: {lt_file}")
                return
                
            # 查找write("Data Atoms")部分的结束位置
            for i in range(atom_section_start, len(lt_content)):
                if 'write(' in line and 'write("Data Atoms")' not in line:
                    atom_section_end = i
                    break
            
            if atom_section_end is None:
                atom_section_end = len(lt_content)
            
            self.logger.info(f"在LT文件中找到原子定义部分: 行 {atom_section_start} 到 {atom_section_end}")
            
            # 收集实际包含原子定义的行
            atom_lines = []
            for i in range(atom_section_start, atom_section_end):
                line = lt_content[i]
                if "$atom:" in line and "." in line:  # 原子定义行的特征
                    atom_lines.append(i)
            
            if len(atom_lines) != len(charges):
                self.logger.warning(f"原子数量({len(atom_lines)})与电荷值数量({len(charges)})不匹配")
                # 如果不匹配，使用较小的数量
                max_update = min(len(atom_lines), len(charges))
                self.logger.info(f"将更新前 {max_update} 个原子的电荷")
            else:
                max_update = len(charges)
            
            # 更新原子电荷
            for i in range(max_update):
                line_idx = atom_lines[i]
                line = lt_content[line_idx]
                
                # 分解原子定义行，格式通常是：
                # $atom:id @atom:type $mol:molid x y z charge
                parts = line.split()
                if len(parts) >= 7:  # 确保行有足够的部分
                    # 最后一部分是电荷
                    parts[-1] = str(charges[i])
                    lt_content[line_idx] = ' '.join(parts) + '\n'
            
            # 写回文件
            with open(lt_file, 'w') as f:
                f.writelines(lt_content)
                
            self.logger.info(f"成功更新LT文件中的电荷值: {lt_file}")
            
        except Exception as e:
            self.logger.error(f"更新LT文件中的电荷时出错: {str(e)}")
            # 记录错误但不抛出异常，使流程可以继续

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
            original_dir = self.work_dir
            
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
            if "original_dir" in locals():
                try:
                    os.chdir(original_dir)
                    self.logger.info(f"恢复原始目录: {original_dir}")
                except Exception as e:
                    self.logger.warning(f"恢复目录失败: {str(e)}")
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
        
        # 动态计算输出频率
        # 目标是在每个阶段获取约100-200帧数据，并考虑到最小和最大输出频率限制

        # 热力学数据输出（保存更多点，因为文件小）
        thermo_freq = max(100, min(equilibration_steps // 200, 5000))

        # NPT阶段轨迹输出（平衡阶段，适当采样）
        trj_freq_npt = max(500, min(equilibration_steps // 100, 10000))

        # NVT阶段轨迹输出（生产阶段，需要足够间隔但不能太稀疏）
        trj_freq_nvt = max(1000, min(production_steps // 150, 20000))

        # RDF输出频率
        rdf_output_freq = max(1000, production_steps // 1000)

        # 打印日志输出便于用户了解
        self.logger.info(f"动态计算的输出频率 - thermo: {thermo_freq}, NPT轨迹: {trj_freq_npt}, NVT轨迹: {trj_freq_nvt}, RDF: {rdf_output_freq}")
        
        # 确定系统名称
        system_name = os.path.basename(system_lt_file).replace('.lt', '')
        
        # 获取系统成分信息 - 从formulation_data获取而不是config
        cations = self.formulation_data.get('cations', [])
        anions = self.formulation_data.get('anions', [])
        solvents = self.formulation_data.get('solvents', [])
        
        # 记录系统成分
        self.logger.info(f"电解液成分 - 阳离子: {len(cations)}个, 阴离子: {len(anions)}个, 溶剂: {len(solvents)}个")
        if cations:
            cation_details = ', '.join([c.get('name', 'Unknown') + '(' + str(c.get('number', 0)) + ')' for c in cations])
            self.logger.info(f"阳离子详情: {cation_details}")
        if anions:
            anion_details = ', '.join([a.get('name', 'Unknown') + '(' + str(a.get('number', 0)) + ')' for a in anions])
            self.logger.info(f"阴离子详情: {anion_details}")
            
        # 从数据文件中读取元素信息
        data_file_path = os.path.join(output_dir, f"{system_name}.data")
        lmpListModifyArr = []
        
        # 定义元素类型函数
        def elementType(element):
            """根据原子质量或注释信息推断元素符号"""
            # 常见元素的原子质量（取整）映射
            mass_map = {
                1: 'H', 4: 'He', 7: 'Li', 9: 'Be', 11: 'B', 12: 'C', 14: 'N', 16: 'O', 19: 'F', 20: 'Ne',
                23: 'Na', 24: 'Mg', 27: 'Al', 28: 'Si', 31: 'P', 32: 'S', 35: 'Cl', 39: 'K', 40: 'Ca',
                45: 'Sc', 48: 'Ti', 51: 'V', 52: 'Cr', 55: 'Mn', 56: 'Fe', 59: 'Co', 58: 'Ni', 63: 'Cu', 
                65: 'Zn', 70: 'Ga', 73: 'Ge', 75: 'As', 79: 'Br', 85: 'Rb', 88: 'Sr', 89: 'Y', 91: 'Zr',
                93: 'Nb', 98: 'Mo', 101: 'Ru', 103: 'Rh', 106: 'Pd', 108: 'Ag', 112: 'Cd', 115: 'In',
                119: 'Sn', 122: 'Sb', 127: 'I', 133: 'Cs', 137: 'Ba', 139: 'La'
            }
            # 尝试从质量映射中查找
            return mass_map.get(element, f"Unknown({element})")
        
        # 读取数据文件，提取元素类型
        atom_type_info = {}  # 存储原子类型、元素和所属分子的映射
        solvent_molecules = []  # 所有溶剂分子列表
        
        try:
            # 先从配方数据获取溶剂分子列表，确保使用前端提交的分子类型
            if self.formulation_data and 'solvents' in self.formulation_data:
                for solvent in self.formulation_data.get('solvents', []):
                    if 'name' in solvent:
                        solvent_name = solvent['name'].replace('+', '').replace('-', '')
                        if solvent_name not in solvent_molecules:
                            solvent_molecules.append(solvent_name)
                
                self.logger.info(f"从formulation_data获取溶剂分子列表: {solvent_molecules}")
            
            # 首先从system.in.list_salt文件获取离子信息
            list_salt_path = os.path.join(output_dir, "system.in.list_salt")
            if os.path.exists(list_salt_path):
                self.logger.info(f"解析离子文件: {list_salt_path}")
                
                # 使用正则表达式匹配行
                group_pattern = re.compile(r'group\s+(\w+)\s+type\s+(.*)')
                variable_pattern = re.compile(r'variable\s+(\w+)_list\s+index\s+"([^"]*)"')
                
                # 存储离子类型和元素
                group_types = {}  # {分子名: [类型列表]}
                group_elements = {}  # {分子名: [元素列表]}
                
                with open(list_salt_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        # 匹配离子组
                        group_match = group_pattern.match(line)
                        if group_match:
                            molecule_name = group_match.group(1)
                            type_str = group_match.group(2)
                            type_list = [int(t) for t in type_str.split()]
                            group_types[molecule_name] = type_list
                            self.logger.info(f"找到离子组 {molecule_name}: 原子类型 {type_list}")
                            continue
                        
                        # 匹配元素列表
                        variable_match = variable_pattern.match(line)
                        if variable_match:
                            mol_name = variable_match.group(1)
                            elements = variable_match.group(2).split()
                            group_elements[mol_name] = elements
                            self.logger.info(f"找到元素列表 {mol_name}: {elements}")
                            continue
                
                # 将原子类型与元素对应
                for mol_name, types in group_types.items():
                    elements = group_elements.get(mol_name, [])
                    
                    for i, atom_type in enumerate(types):
                        # 确定元素，如果元素不够，循环使用
                        element = elements[i % len(elements)] if elements else "Unknown"
                        
                        # 添加到映射
                        atom_type_info[atom_type] = {
                            'element': element,
                            'molecule': mol_name
                        }
                        self.logger.info(f"映射离子: 类型 {atom_type} -> 元素 {element} 在分子 {mol_name}")
            
            # 然后从system.data文件读取所有原子类型
            with open(data_file_path, 'r', encoding='utf-8') as f:
                content = f.read()
                
                # 查找Masses段落
                masses_match = re.search(r'Masses\s*\n\s*\n(.*?)(?=\n\s*\n\w+|\Z)', content, re.DOTALL)
                if masses_match:
                    masses_section = masses_match.group(1)
                    self.logger.info(f"找到Masses段落，长度: {len(masses_section)} 字符")
                    
                    # 解析每一行
                    type1_positions = []  # 记录type1出现的位置
                    
                    for line in masses_section.strip().split('\n'):
                        line = line.strip()
                        if not line:
                            continue
                            
                        # 分离数据和注释
                        parts = line.split('#', 1)
                        if len(parts) < 2:
                            continue
                            
                        data = parts[0].strip()
                        comment = parts[1].strip()
                        
                        # 解析数据部分
                        items = data.split()
                        if len(items) < 2:
                            continue
                            
                        atom_type = int(items[0])
                        mass = float(items[1])
                        atomic_mass = round(mass)
                        
                        # 确定元素:
                        # 1. 首先检查是否已经在离子映射中
                        if atom_type in atom_type_info:
                            element = atom_type_info[atom_type]['element']
                        # 2. 否则，检查注释中是否有元素名称 (非typeX格式)
                        elif not comment.startswith('type'):
                            # 注释第一个单词可能是元素符号
                            element = comment.split()[0]
                        # 3. 最后，尝试从质量推断
                        else:
                            element = elementType(atomic_mass)
                        
                        # 记录原子类型和元素
                        if atom_type not in atom_type_info:
                            atom_type_info[atom_type] = {
                                'mass': mass,
                                'element': element,
                                'comment': comment,
                                'molecule': None  # 溶剂分子将在后续步骤中分配
                            }
                        
                        # 检测type1标记 (用于识别溶剂开始位置)
                        if 'type1' in comment and atom_type > 9:  # 排除离子类型
                            type1_positions.append(atom_type)
                            self.logger.info(f"在原子类型 {atom_type} 处检测到溶剂标记 (type1)")
                    
                    # 识别溶剂 - 如果从formulation_data没有获取到溶剂信息，则从文件检测
                    if not solvent_molecules:
                        self.logger.info("未从formulation_data获取到溶剂信息，尝试从文件检测")
                        for solvent in ['EC', 'PC', 'DMC', 'EMC', 'DEC', 'ACN', 'DMSO']:
                            if os.path.exists(os.path.join(output_dir, f"{solvent}.lt")):
                                solvent_molecules.append(solvent)
                        
                        self.logger.info(f"从文件检测到的溶剂分子: {solvent_molecules}")
                    
                    # 根据type1标记分配溶剂类型
                    if type1_positions and solvent_molecules:
                        self.logger.info(f"找到 {len(type1_positions)} 个type1标记: {type1_positions}")
                        
                        # 处理溶剂映射
                        if len(solvent_molecules) > 0:
                            current_solvent_index = 0
                            current_solvent = solvent_molecules[current_solvent_index]
                            molecule_types = {current_solvent: []}
                            
                            # 查找连续的typeX标记，作为一个溶剂分子
                            in_solvent = False
                            curr_type_num = 0
                            prev_type_num = 0
                            
                            # 排序原子类型，确保按顺序处理
                            sorted_types = sorted(atom_type_info.keys())
                            for atom_type in sorted_types:
                                info = atom_type_info[atom_type]
                                if 'comment' not in info:
                                    continue
                                    
                                comment = info['comment']
                                
                                # 检查是否是离子类型，已经有分子归属的跳过
                                if info.get('molecule') and info.get('molecule') in group_types.keys():
                                    continue
                                
                                # 检查是否包含typeX标记
                                type_match = re.search(r'type(\d+)', comment)
                                if type_match:
                                    curr_type_num = int(type_match.group(1))
                                    
                                    # 检测新溶剂的开始 (type1)
                                    if curr_type_num == 1 and not in_solvent:
                                        in_solvent = True
                                        if current_solvent_index < len(solvent_molecules):
                                            current_solvent = solvent_molecules[current_solvent_index]
                                            molecule_types[current_solvent] = []
                                        else:
                                            # 如果溶剂已经用完，使用Unknown
                                            current_solvent = f"Unknown_{current_solvent_index}"
                                            molecule_types[current_solvent] = []
                                    
                                    # 检测溶剂分子的结束 (如果当前类型比前一个小，说明是新分子)
                                    if curr_type_num < prev_type_num and curr_type_num == 1:
                                        # 结束当前溶剂，开始新溶剂
                                        current_solvent_index += 1
                                        if current_solvent_index < len(solvent_molecules):
                                            current_solvent = solvent_molecules[current_solvent_index]
                                            molecule_types[current_solvent] = []
                                        else:
                                            # 如果溶剂已经用完，使用Unknown
                                            current_solvent = f"Unknown_{current_solvent_index}"
                                            molecule_types[current_solvent] = []
                                    
                                    # 将当前原子类型添加到当前溶剂
                                    if in_solvent:
                                        molecule_types[current_solvent].append(atom_type)
                                        atom_type_info[atom_type]['molecule'] = current_solvent
                                    
                                    prev_type_num = curr_type_num
                            
                            # 输出溶剂分配情况
                            self.logger.info("\n溶剂原子类型分配:")
                            for solvent, types in molecule_types.items():
                                self.logger.info(f"溶剂 {solvent}: {len(types)} 个原子类型 - {types}")
        except Exception as e:
            self.logger.error(f"读取数据文件失败: {str(e)}")
            import traceback
            traceback.print_exc()
            # 如果读取失败，使用默认元素列表
            atom_type_info = {}
            for i, element in enumerate(["Li", "P", "F", "F", "F", "F", "F", "F", "C", "C", "C", "O", "O", "H", "H", "H", "H"], start=1):
                atom_type_info[i] = {'element': element, 'molecule': 'Default'}
            self.logger.warning(f"使用默认元素列表")
        
        # 准备element_list
        lmpListModifyArr = []
        for atom_type in sorted(atom_type_info.keys()):
            element = atom_type_info[atom_type]['element']
            lmpListModifyArr.append(element)
        
        # 生成元素列表文件
        element_list = " ".join(lmpListModifyArr)
        self.logger.info(f"元素列表: {element_list}")
        
        # 写入元素列表到.in.list文件
        list_file_path = os.path.join(output_dir, f"{system_name}.in.list")
        with open(list_file_path, 'w', encoding='utf-8') as f:
            f.write(f"variable element_list index \"{element_list}\"\n")
        
        # 创建RDF对的信息
        findLiO = ''
        fixRdf = ''
        
        # 定义阳离子和阴离子元素
        cation_elements = ['Li', 'Na', 'K', 'Mg', 'Ca', 'Zn', 'Al']
        anion_elements = ['O', 'N', 'P', 'B', 'F', 'S', 'Cl']
        
        self.logger.info(f"开始生成RDF对，溶剂分子列表: {solvent_molecules}")
        self.logger.info(f"可用的原子类型信息: {len(atom_type_info)}个类型")
        
        # 输出原子类型信息摘要
        atom_info_summary = []
        for atom_type, info in sorted(atom_type_info.items())[:10]:  # 只显示前10个，避免日志过长
            molecule = info.get('molecule', 'Unknown')
            element = info.get('element', 'Unknown')
            atom_info_summary.append(f"类型{atom_type}({element}@{molecule})")
        
        self.logger.info(f"原子类型信息示例: {', '.join(atom_info_summary)}...")
        
        # 遍历原子类型信息创建RDF对
        for atom_type, info in sorted(atom_type_info.items()):
            element = info['element']
            molecule = info.get('molecule', 'Unknown')
            
            if element in cation_elements:
                for anion_type, anion_info in sorted(atom_type_info.items()):
                    anion_element = anion_info['element']
                    anion_molecule = anion_info.get('molecule', 'Unknown')
                    
                    if anion_element in anion_elements:
                        # 对于离子-离子对，使用 rdf_cation_anion 格式
                        if molecule in group_types.keys() and anion_molecule in group_types.keys():
                            findLiO += f'{atom_type} {anion_type} '
                            rdf_label = f'rdf_{element}_{anion_element}_{anion_molecule}'
                            cn_label = f'cn_{element}_{anion_element}_{anion_molecule}'
                            fixRdf += f'{rdf_label} {cn_label} '
                            self.logger.info(f"添加离子-离子RDF对: {rdf_label} (类型: {atom_type}-{anion_type})")
                        # 对于离子-溶剂对，使用 rdf_cation_anion_solventname 格式
                        elif molecule in group_types.keys() and anion_molecule not in group_types.keys() and anion_molecule in solvent_molecules:
                            findLiO += f'{atom_type} {anion_type} '
                            rdf_label = f'rdf_{element}_{anion_element}_{anion_molecule}'
                            cn_label = f'cn_{element}_{anion_element}_{anion_molecule}'
                            fixRdf += f'{rdf_label} {cn_label} '  
                            self.logger.info(f"添加离子-溶剂RDF对: {rdf_label} (类型: {atom_type}-{anion_type}, 溶剂: {anion_molecule})")
        
        self.logger.info(f"RDF对生成完成，共 {len(findLiO.split()) // 2} 个对")
        
        # 创建LAMMPS输入文件内容
        in_file_content = f"""# LAMMPS输入脚本 - 生成于 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
# ----------------- Variable Section -----------------
variable infile string {system_name}
variable outname string {system_name}
include {system_name}.in.init
read_data {system_name}.data
include {system_name}.in.settings

# 元素列表
variable element_list index "{element_list}"
variable rdf_pair string "{findLiO}"

# 模拟参数
variable Nsteps_NPT equal {equilibration_steps}
variable Nsteps_NVT equal {production_steps}
variable Freq_trj_npt equal {trj_freq_npt}
variable Freq_trj_nvt equal {trj_freq_nvt}
variable thermo_freq equal {thermo_freq}
variable timestep equal {time_step}
variable Temp_NPT equal {temp}
variable Temp_NVT equal {temp}
variable rdf_output_freq equal {rdf_output_freq}
variable pressure equal {pressure}
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
fix rdff1 all ave/time ${{rdf_output_freq}} 1000 ${{Nsteps_NVT}} c_rdfc1[*] file out_rdf.dat mode vector title3 "RDF {fixRdf}"

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
export PATH=/public/software/lammps/mpich_install/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/software/lammps/mpich_install/lib:/public/software/intel/oneapi/compiler/2024.0/lib
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
        
        # 对molecules_info进行去重处理，避免同一分子被添加多次
        deduplicated_molecules = {}
        for molecule in molecules_info:
            name = molecule.get("name", "")
            if not name:
                continue
                
            # 如果这个分子已经存在，合并数量
            if name in deduplicated_molecules:
                existing_count = deduplicated_molecules[name].get("number", deduplicated_molecules[name].get("count", 0))
                new_count = molecule.get("number", molecule.get("count", 0))
                
                # 合并数量
                deduplicated_molecules[name]["number"] = existing_count + new_count
                self.logger.warning(f"检测到重复分子 {name}，合并数量: {existing_count} + {new_count} = {existing_count + new_count}")
            else:
                # 新分子
                deduplicated_molecules[name] = molecule
        
        # 将字典转回列表
        molecules_info = list(deduplicated_molecules.values())
        self.logger.info(f"去重后的分子列表包含 {len(molecules_info)} 个分子")
        
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
        """使用RESP2计算分子的电荷，在任务特定的临时目录中执行以避免多任务冲突"""
        
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
            
            # 创建任务特定的临时RESP目录
            task_resp_dir = os.path.join(self.work_dir, "temp_resp", mol_name)
            os.makedirs(task_resp_dir, exist_ok=True)
            self.logger.info(f"创建任务特定的临时RESP目录: {task_resp_dir}")
            
            # 确保电荷保存目录存在
            charge_save_dir = self.initial_solvent_path
            os.makedirs(charge_save_dir, exist_ok=True)
            
            # 复制PDB文件到临时RESP目录
            dest_pdb = os.path.join(task_resp_dir, f"{mol_name}.pdb")
            shutil.copy2(pdb_file, dest_pdb)
            self.logger.info(f"已复制PDB文件到临时RESP目录: {pdb_file} -> {dest_pdb}")
            
            # 从全局RESP目录复制必要的脚本文件到临时目录
            resp_script_source = os.path.join(resp_path, "RESP2.sh")
            if not os.path.exists(resp_script_source):
                self.logger.error(f"RESP2脚本不存在于全局目录: {resp_script_source}")
                return
                
            # 复制RESP2.sh和其他必要文件到临时目录
            for script_file in ["RESP2.sh", "gen_esp_dat.py", "resp2.py"]:
                src_file = os.path.join(resp_path, script_file)
                if os.path.exists(src_file):
                    dest_file = os.path.join(task_resp_dir, script_file)
                    shutil.copy2(src_file, dest_file)
                    self.logger.info(f"已复制RESP脚本文件: {src_file} -> {dest_file}")
            
            # 在临时目录中运行RESP脚本
            resp_script = os.path.join(task_resp_dir, "RESP2.sh")
            
            # 切换到临时RESP目录
            original_dir = os.getcwd()
            os.chdir(task_resp_dir)
            self.logger.info(f"切换到临时RESP目录: {task_resp_dir}")
            
            # 修改RESP2.sh文件，设置Gaussian绝对路径
            g16_path = "/public/software/g16/"
            
            # 读取RESP2.sh内容
            with open(resp_script, 'r') as f:
                script_content = f.read()
            
            # 替换Gaussian路径
            script_content = script_content.replace('Gaussian=g16', f'Gaussian="{g16_path}/g16"')
            script_content = script_content.replace('Gaussian=g09', f'Gaussian="{g16_path}/g16"')
            
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
            self.logger.info(f"在临时目录中执行RESP命令: {cmd}")
            
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            stdout, stderr = process.communicate()
            
            # 检查命令执行是否成功
            if process.returncode != 0:
                self.logger.error(f"RESP命令执行失败: {stderr}")
                self.logger.debug(f"RESP stdout: {stdout}")
                # 恢复原始目录
                os.chdir(original_dir)
                return
            
            self.logger.info(f"RESP命令执行成功")
            
            # 查找临时目录中生成的电荷文件
            charge_file = os.path.join(task_resp_dir, f"{mol_name}.chg")
            if os.path.exists(charge_file):
                # 复制电荷文件到initial_solvent目录
                dest_charge = os.path.join(charge_save_dir, f"{mol_name}.chg")
                shutil.copy2(charge_file, dest_charge)
                self.logger.info(f"已将电荷文件复制到initial_solvent目录: {charge_file} -> {dest_charge}")
                
                # 立即更新existing_molecules.json
                self._update_exist_solvent_json(mol_name)
            else:
                self.logger.warning(f"RESP未生成电荷文件: {charge_file}")
                # 搜索可能使用其他名称生成的电荷文件
                chg_files = [f for f in os.listdir(task_resp_dir) if f.endswith('.chg')]
                if chg_files:
                    self.logger.info(f"找到其他电荷文件: {chg_files}")
                    # 使用第一个找到的电荷文件
                    first_chg = os.path.join(task_resp_dir, chg_files[0])
                    dest_charge = os.path.join(charge_save_dir, f"{mol_name}.chg")
                    shutil.copy2(first_chg, dest_charge)
                    self.logger.info(f"已复制找到的电荷文件: {first_chg} -> {dest_charge}")
            
            # 恢复原始目录
            os.chdir(original_dir)
            self.logger.info(f"恢复原始工作目录: {original_dir}")
            
            # 清理临时目录中的大文件以节省空间，但保留关键文件以便调试
            try:
                for root, dirs, files in os.walk(task_resp_dir):
                    for file in files:
                        if file.endswith(('.fchk', '.chk', '.log')) and not file.startswith(mol_name):
                            file_path = os.path.join(root, file)
                            if os.path.getsize(file_path) > 10*1024*1024:  # 大于10MB
                                os.remove(file_path)
                                self.logger.info(f"删除大型临时文件: {file_path}")
            except Exception as e:
                self.logger.warning(f"清理临时文件时出错: {str(e)}")
            
            self.logger.info(f"RESP电荷计算完成: {mol_name}")
        
        except Exception as e:
            self.logger.error(f"生成RESP电荷时出错: {str(e)}")
            # 记录更详细的异常信息
            import traceback
            self.logger.error(f"详细错误信息: {traceback.format_exc()}")

    def _update_exist_solvent_json(self, mol_name: str) -> None:
        """更新existing_molecules.json文件，保存溶剂的SMILE表达式便于后续查找
        
        Args:
            mol_name: 分子名称
        """
        try:
            # 获取分子的SMILE字符串
            smile = ""
            for solvent in self.formulation_data.get("solvents", []):
                if solvent.get("name") == mol_name:
                    smile = solvent.get("smile", "")
                    break
            
            if not smile:
                self.logger.warning(f"未找到分子 {mol_name} 的SMILE字符串")
                return
            
            # 读取现有的existing_molecules.json
            existing_molecules_file = os.path.join(self.initial_solvent_path, "existing_molecules.json")
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
            
            # 设置溶剂文件路径
            pdb_file_path = os.path.join(self.initial_solvent_path, f"{mol_name}.pdb")
            lt_file_path = os.path.join(self.initial_solvent_path, f"{mol_name}.lt")
            chg_file_path = os.path.join(self.initial_solvent_path, f"{mol_name}.chg")
            
            # 更新数据 - 添加到溶剂部分
            existing_data["solvent"][mol_name] = {
                "name": mol_name,
                "smile": smile,
                "pdb_file": pdb_file_path if os.path.exists(pdb_file_path) else "",
                "lt_file": lt_file_path if os.path.exists(lt_file_path) else "",
                "charge_file": chg_file_path if os.path.exists(chg_file_path) else "",
                "last_updated": datetime.now().isoformat()
            }
            
            # 保存更新后的数据
            with open(existing_molecules_file, 'w') as f:
                json.dump(existing_data, f, indent=2)
            
            self.logger.info(f"已更新existing_molecules.json，添加溶剂分子: {mol_name}，SMILE: {smile}")
            
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

    def analyze_trajectory(self, trajectory_file=None, log_file=None):
        """
        分析模拟轨迹文件，计算物理性质，并保存结果
        
        Args:
            trajectory_file: 轨迹文件路径，默认为None
            log_file: 日志文件路径，默认为None
        
        Returns:
            dict: 分析结果字典
        """
        try:
            import os
            import sys
            import json
            import subprocess
            from pathlib import Path
            
            # 获取分析脚本的绝对路径
            aisuan_root = Path(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            analysis_script = str(aisuan_root / "analyze_electrolyte_results.py")
            
            if not os.path.exists(analysis_script):
                self.logger.error(f"分析脚本不存在: {analysis_script}")
                return {"status": "failed", "error": f"分析脚本不存在: {analysis_script}"}
            
            # 确保脚本可执行
            os.chmod(analysis_script, 0o755)
            
            # 准备分析参数
            task_dir = os.path.dirname(self.output_dir) if self.output_dir else self.work_dir
            output_dir = os.path.join(task_dir, "analysis")
            
            # 从配方数据中获取温度和浓度
            temperature = 298.15  # 默认温度 (K)
            concentration = 1.0   # 默认浓度 (mol/L)
            
            if hasattr(self, 'formulation_data') and self.formulation_data:
                if 'parameters' in self.formulation_data:
                    temperature = self.formulation_data['parameters'].get('temperature', temperature)
                    
                # 尝试从盐的摩尔浓度计算总浓度
                if 'salts' in self.formulation_data:
                    total_salt_concentration = 0
                    for salt in self.formulation_data['salts']:
                        if 'concentration' in salt:
                            total_salt_concentration += salt['concentration']
                    
                    if total_salt_concentration > 0:
                        concentration = total_salt_concentration
            
            self.logger.info(f"使用参数进行分析: 温度={temperature}K, 浓度={concentration}mol/L")
            
            # 构建命令
            cmd = [
                "python", analysis_script,
                "--task_dir", task_dir,
                "--output_dir", output_dir,
                "--analysis_type", "all",
                "--temperature", str(temperature),
                "--concentration", str(concentration)
            ]
            
            # 执行分析脚本
            self.logger.info(f"执行分析命令: {' '.join(cmd)}")
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                self.logger.error(f"分析脚本执行失败: {stderr}")
                return {"status": "failed", "error": stderr}
            
            self.logger.info(f"分析脚本执行成功")
            self.logger.debug(f"分析脚本输出: {stdout}")
            
            # 读取分析结果
            results = {"status": "success"}
            
            # 读取MSD分析结果
            msd_results_file = os.path.join(output_dir, "msd", "msd_analysis_results.txt")
            if os.path.exists(msd_results_file):
                with open(msd_results_file, 'r') as f:
                    msd_text = f.read()
                    
                    # 解析扩散系数
                    diffusion_line = next((line for line in msd_text.split('\n') if '扩散系数' in line), None)
                    if diffusion_line:
                        try:
                            diffusion_value = float(diffusion_line.split(':')[1].strip().split()[0])
                            results["diffusion_coefficient"] = diffusion_value
                        except:
                            self.logger.warning("无法解析扩散系数")
                    
                    # 解析离子电导率
                    conductivity_line = next((line for line in msd_text.split('\n') if '离子电导率' in line), None)
                    if conductivity_line:
                        try:
                            conductivity_value = float(conductivity_line.split(':')[1].strip().split()[0])
                            results["ionic_conductivity"] = conductivity_value
                        except:
                            self.logger.warning("无法解析离子电导率")
            
            # 获取图表路径
            msd_plot = os.path.join(output_dir, "msd", "msd_plot.png")
            rdf_plot = os.path.join(output_dir, "rdf", "rdf_plot.png")
            
            if os.path.exists(msd_plot):
                results["plots"] = results.get("plots", []) + [{"name": "均方位移 (MSD)", "path": msd_plot}]
            
            if os.path.exists(rdf_plot):
                results["plots"] = results.get("plots", []) + [{"name": "径向分布函数 (RDF)", "path": rdf_plot}]
            
            # 整合结果
            results["results"] = {
                "diffusion_coefficient": results.get("diffusion_coefficient"),
                "ionic_conductivity": results.get("ionic_conductivity"),
                "temperature": temperature,
                "concentration": concentration,
                "analysis_path": output_dir
            }
            
            self.logger.info(f"分析结果: {json.dumps(results, indent=2)}")
            return results
        
        except Exception as e:
            self.logger.error(f"分析轨迹时出错: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            return {"status": "failed", "error": str(e)}

    def check_calculation_completion(self):
        """
        检查LAMMPS计算是否完成，如果完成则自动运行分析
        
        Returns:
            bool: 如果计算完成返回True，否则返回False
        """
        import os
        import time
        
        if not hasattr(self, 'output_dir') or not self.output_dir:
            self.logger.warning("没有设置输出目录，无法检查计算完成状态")
            return False
        
        # 检查是否存在输出目录
        if not os.path.exists(self.output_dir):
            self.logger.debug(f"输出目录不存在: {self.output_dir}")
            return False
        
        # 检查是否存在标准LAMMPS输出文件
        log_file = os.path.join(self.output_dir, "log.lammps")
        if not os.path.exists(log_file):
            self.logger.debug(f"LAMMPS日志文件不存在: {log_file}")
            return False
        
        # 检查轨迹文件
        trajectory_file = os.path.join(self.output_dir, "NVT_system.lammpstrj")
        if not os.path.exists(trajectory_file):
            self.logger.debug(f"轨迹文件不存在: {trajectory_file}")
            return False
        
        # 检查log文件中是否包含计算完成的标记
        with open(log_file, 'r') as f:
            log_content = f.read()
            if "Total wall time:" not in log_content:
                self.logger.debug("LAMMPS计算尚未完成")
                return False
        
        self.logger.info("LAMMPS计算已完成")
        
        # 如果已经存在分析结果，不再重复分析
        analysis_dir = os.path.join(os.path.dirname(self.output_dir), "analysis")
        if os.path.exists(analysis_dir):
            self.logger.info(f"分析目录已存在: {analysis_dir}，跳过自动分析")
            return True
        
        # 自动运行分析
        self.logger.info("开始自动分析轨迹数据")
        analysis_result = self.analyze_trajectory(
            trajectory_file=trajectory_file,
            log_file=log_file
        )
        
        if analysis_result.get("status") == "success":
            self.logger.info("自动分析完成")
        else:
            self.logger.warning(f"自动分析失败: {analysis_result.get('error', '未知错误')}")
        
        return True

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