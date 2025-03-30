"""
分子库管理模块

该模块提供分子库的管理功能，包括获取分子文件、添加新分子等。
"""

import os
import json
import shutil
import logging
from typing import Dict, List, Optional, Union, Any

# 设置日志
logger = logging.getLogger("MoleculeLibrary")

class MoleculeLibrary:
    """分子库管理类，用于管理和访问分子文件"""
    
    def __init__(self, library_path: str):
        """初始化分子库
        
        Args:
            library_path: 分子库根目录路径
        """
        self.library_path = library_path
        self.metadata_path = os.path.join(library_path, "metadata")
        self.ions_path = os.path.join(library_path, "ions")
        self.solvents_path = os.path.join(library_path, "solvents")
        
        # 加载元数据
        self.molecules = self._load_json(os.path.join(self.metadata_path, "molecules.json"))
        self.smiles_index = self._load_json(os.path.join(self.metadata_path, "smiles_index.json"))
        
        # 设置日志
        self._setup_logger()
    
    def _setup_logger(self):
        """设置日志"""
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    def _load_json(self, file_path: str) -> Dict:
        """加载JSON文件
        
        Args:
            file_path: JSON文件路径
            
        Returns:
            Dict: 解析后的JSON内容，如果文件不存在则返回空字典
        """
        if not os.path.exists(file_path):
            logger.warning(f"文件不存在: {file_path}")
            return {}
        
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            logger.error(f"JSON解析错误: {file_path}")
            return {}
        except Exception as e:
            logger.error(f"加载文件时出错: {file_path}, 错误: {str(e)}")
            return {}
    
    def _save_json(self, data: Dict, file_path: str) -> bool:
        """保存JSON文件
        
        Args:
            data: 要保存的数据
            file_path: 目标文件路径
            
        Returns:
            bool: 是否成功保存
        """
        try:
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            return True
        except Exception as e:
            logger.error(f"保存文件时出错: {file_path}, 错误: {str(e)}")
            return False
    
    def _get_molecule_path(self, name: str, is_cation: bool = False, is_anion: bool = False, is_solvent: bool = False) -> str:
        """获取分子文件目录路径
        
        Args:
            name: 分子名称
            is_cation: 是否是阳离子
            is_anion: 是否是阴离子
            is_solvent: 是否是溶剂
                    
        Returns:
            str: 分子文件目录路径
        """
        # 如果未指定类别标志，则自动判断
        if not any([is_cation, is_anion, is_solvent]):
            if name in self.molecules.get("ions", {}).get("cations", []):
                is_cation = True
            elif name in self.molecules.get("ions", {}).get("anions", []):
                is_anion = True
            elif name in self.molecules.get("solvents", []):
                is_solvent = True
            else:
                logger.warning(f"未知分子: {name}")
                return ""
        
        # 根据类别返回路径
        if is_cation:
            return os.path.join(self.ions_path, "cations", name)
        elif is_anion:
            return os.path.join(self.ions_path, "anions", name)
        elif is_solvent:
            return os.path.join(self.solvents_path, name)
        else:
            logger.warning(f"未指定有效的分子类别")
            return ""
    
    def get_molecule(self, name: str, category: str = None, target_dir: str = None, use_symlink: bool = True) -> Dict[str, str]:
        """获取分子文件
        
        Args:
            name: 分子名称
            category: 分子类别，可以是'cations', 'anions'或'solvents'
            target_dir: 目标目录，如果为None则不复制
            use_symlink: 是否使用符号链接，默认为True
            
        Returns:
            Dict: 包含文件路径的字典，如果失败则返回空字典
        """
        # 判断分子类别
        is_cation = False
        is_anion = False
        is_solvent = False
        
        if category:
            is_cation = category == "cations"
            is_anion = category == "anions"
            is_solvent = category == "solvents"
        else:
            # 尝试在所有类别中查找
            if name in self.molecules.get("ions", {}).get("cations", []):
                is_cation = True
            elif name in self.molecules.get("ions", {}).get("anions", []):
                is_anion = True
            elif name in self.molecules.get("solvents", []):
                is_solvent = True
        
        # 获取分子目录
        mol_dir = self._get_molecule_path(name, is_cation, is_anion, is_solvent)
        if not mol_dir:
            return {}
            
        # 获取文件路径
        result_files = {}
        
        # 离子需要pdb和lt文件
        if is_cation or is_anion:
            for ext in ["pdb", "lt"]:
                src_file = os.path.join(mol_dir, f"{name}.{ext}")
                if os.path.exists(src_file):
                    result_files[ext] = src_file
                else:
                    logger.warning(f"离子 {name} 的 {ext} 文件不存在: {src_file}")
        
        # 溶剂需要pdb、lt和chg文件
        elif is_solvent:
            for ext in ["pdb", "lt", "chg"]:
                src_file = os.path.join(mol_dir, f"{name}.{ext}")
                if os.path.exists(src_file):
                    result_files[ext] = src_file
                else:
                    logger.warning(f"溶剂 {name} 的 {ext} 文件不存在: {src_file}")
        
        # 如果没有指定目标目录，直接返回文件路径
        if target_dir is None:
            return result_files
        
        # 确保目标目录存在
        os.makedirs(target_dir, exist_ok=True)
        
        # 复制或链接文件到目标目录
        target_files = {}
        
        # 检查必要的文件是否存在
        required_exts = ["pdb", "lt"]
        if is_solvent:
            required_exts.append("chg")
            
        missing_files = [ext for ext in required_exts if ext not in result_files]
        if missing_files:
            logger.error(f"分子 {name} 缺少必要的文件: {', '.join(missing_files)}")
            # 尽管缺少文件，仍然尝试处理可用文件
        
        for ext, src_file in result_files.items():
            target_file = os.path.join(target_dir, f"{name}.{ext}")
            
            # 如果文件已存在，先删除
            if os.path.exists(target_file):
                try:
                    os.remove(target_file)
                except Exception as e:
                    logger.warning(f"删除已存在的文件时出错: {target_file}, 错误: {str(e)}")
                    continue
            
            try:
                # 验证源文件可读性
                try:
                    with open(src_file, 'rb') as test_file:
                        pass  # 仅测试文件可读性
                except Exception as e:
                    logger.error(f"源文件不可读: {src_file}, 错误: {str(e)}")
                    raise IOError(f"源文件不可读: {src_file}")
                
                if use_symlink:
                    # 使用符号链接
                    os.symlink(src_file, target_file)
                    logger.info(f"已创建符号链接: {target_file} -> {src_file}")
                else:
                    # 复制文件
                    shutil.copy2(src_file, target_file)
                    logger.info(f"已复制文件: {src_file} -> {target_file}")
                
                # 验证目标文件存在和可读性
                if not os.path.exists(target_file):
                    logger.error(f"目标文件未创建成功: {target_file}")
                    raise FileNotFoundError(f"目标文件未创建成功: {target_file}")
                
                try:
                    with open(target_file, 'rb') as test_file:
                        pass  # 测试文件可读性
                except Exception as e:
                    logger.error(f"目标文件不可读: {target_file}, 错误: {str(e)}")
                    raise IOError(f"目标文件不可读: {target_file}")
                
                target_files[ext] = target_file
            except Exception as e:
                logger.error(f"处理文件时出错: {src_file} -> {target_file}, 错误: {str(e)}")
        
        return target_files
    
    def find_by_smile(self, smile: str) -> str:
        """通过SMILE字符串查找分子
        
        Args:
            smile: SMILE字符串
            
        Returns:
            str: 分子名称，如果未找到则返回空字符串
        """
        return self.smiles_index.get(smile, "")
    
    def add_molecule(self, name: str, files: Dict[str, str], category: str, smile: str = None) -> bool:
        """添加新分子到库中
        
        Args:
            name: 分子名称
            files: 文件路径字典，格式为{ext: path}
            category: 分子类别，可以是'cations', 'anions'或'solvents'
            smile: SMILE字符串，仅对溶剂需要
            
        Returns:
            bool: 是否成功添加
        """
        # 检查类别是否有效
        if category not in ["cations", "anions", "solvents"]:
            logger.error(f"无效的分子类别: {category}")
            return False
        
        # 溶剂必须提供SMILE
        if category == "solvents" and not smile:
            logger.error(f"添加溶剂时必须提供SMILE字符串: {name}")
            return False
        
        # 创建分子目录
        mol_dir = ""
        if category == "cations":
            mol_dir = os.path.join(self.ions_path, "cations", name)
            # 更新元数据
            if name not in self.molecules.get("ions", {}).get("cations", []):
                self.molecules.setdefault("ions", {}).setdefault("cations", []).append(name)
        elif category == "anions":
            mol_dir = os.path.join(self.ions_path, "anions", name)
            # 更新元数据
            if name not in self.molecules.get("ions", {}).get("anions", []):
                self.molecules.setdefault("ions", {}).setdefault("anions", []).append(name)
        elif category == "solvents":
            mol_dir = os.path.join(self.solvents_path, name)
            # 更新元数据
            if name not in self.molecules.get("solvents", []):
                self.molecules.setdefault("solvents", []).append(name)
            # 更新SMILE索引
            self.smiles_index[smile] = name
        
        # 创建目录
        os.makedirs(mol_dir, exist_ok=True)
        
        # 复制文件
        success = True
        for ext, src_file in files.items():
            if not os.path.exists(src_file):
                logger.warning(f"源文件不存在: {src_file}")
                success = False
                continue
            
            target_file = os.path.join(mol_dir, f"{name}.{ext}")
            try:
                shutil.copy2(src_file, target_file)
                logger.info(f"已复制文件: {src_file} -> {target_file}")
                
                # 对于溶剂，添加SMILE文件
                if category == "solvents" and ext == "pdb" and smile:
                    smile_file = os.path.join(mol_dir, f"{name}.smile")
                    with open(smile_file, 'w') as f:
                        f.write(smile)
                    logger.info(f"已创建SMILE文件: {smile_file}")
            except Exception as e:
                logger.error(f"复制文件时出错: {src_file} -> {target_file}, 错误: {str(e)}")
                success = False
        
        # 保存元数据
        if success:
            self._save_json(self.molecules, os.path.join(self.metadata_path, "molecules.json"))
            if category == "solvents":
                self._save_json(self.smiles_index, os.path.join(self.metadata_path, "smiles_index.json"))
        
        return success
    
    def list_molecules(self, category: str = None) -> Dict:
        """列出分子库中的分子
        
        Args:
            category: 可选的过滤类别，可以是'cations', 'anions'或'solvents'
            
        Returns:
            Dict: 分子列表，按类别组织
        """
        if category is None:
            return self.molecules
        elif category == "cations":
            return {"cations": self.molecules.get("ions", {}).get("cations", [])}
        elif category == "anions":
            return {"anions": self.molecules.get("ions", {}).get("anions", [])}
        elif category == "solvents":
            return {"solvents": self.molecules.get("solvents", [])}
        else:
            logger.warning(f"未知类别: {category}")
            return {}
    
    def get_molecule_path(self, molecule_name):
        """
        公共方法：获取分子文件目录路径
        :param molecule_name: 分子名称
        :return: 文件路径
        """
        return self._get_molecule_path(molecule_name) 