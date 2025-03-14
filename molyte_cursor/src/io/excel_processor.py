"""
Excel数据处理器模块

该模块提供了读取和处理Excel数据的功能，包括自动计算分子数量
"""

import os
import logging
import re
from typing import Dict, List, Tuple, Optional, Union, Any
import pandas as pd
from openpyxl import load_workbook

from ..utils.molecule_calculator import MoleculeCalculator

class ExcelProcessor:
    """
    Excel数据处理器
    
    处理Excel输入数据，自动计算分子数量，避免用户手动计算
    """
    
    def __init__(self, logger=None):
        """
        初始化Excel处理器
        
        Args:
            logger: 日志记录器实例
        """
        self.logger = logger or logging.getLogger(__name__)
        self.molecule_calculator = MoleculeCalculator(logger)
    
    def read_excel(self, file_path: str, sheet_name: Optional[str] = None, header_row: int = 1) -> pd.DataFrame:
        """
        读取Excel文件并返回DataFrame
        
        Args:
            file_path: Excel文件路径
            sheet_name: 工作表名称，如果为None则使用第一个工作表
            header_row: 表头行号（从1开始）
            
        Returns:
            包含Excel数据的DataFrame
        """
        self.logger.info(f"读取Excel文件: {file_path}")
        
        try:
            # 使用pandas读取Excel
            df = pd.read_excel(file_path, sheet_name=sheet_name, header=header_row-1)
            self.logger.debug(f"成功读取Excel，包含 {len(df)} 行数据")
            return df
        except Exception as e:
            self.logger.error(f"读取Excel文件出错: {e}")
            raise
    
    def get_headers(self, file_path: str, sheet_name: Optional[str] = None, header_row: int = 1) -> List[str]:
        """
        获取Excel文件的表头
        
        Args:
            file_path: Excel文件路径
            sheet_name: 工作表名称，如果为None则使用第一个工作表
            header_row: 表头行号（从1开始）
            
        Returns:
            表头列表
        """
        self.logger.info(f"获取Excel表头: {file_path}")
        
        try:
            # 使用openpyxl读取表头
            wb = load_workbook(file_path, read_only=True)
            if sheet_name is None:
                sheet = wb.active
            else:
                sheet = wb[sheet_name]
            
            headers = [cell.value for cell in sheet[header_row] if cell.value is not None]
            self.logger.debug(f"获取到表头: {headers}")
            return headers
        except Exception as e:
            self.logger.error(f"获取Excel表头出错: {e}")
            raise
    
    def identify_columns(self, headers: List[str]) -> Tuple[List[str], List[str], List[str], Optional[str]]:
        """
        根据表头识别阳离子、阴离子和溶剂的列
        
        Args:
            headers: 表头列表
            
        Returns:
            阳离子列、阴离子列、溶剂列和名称列的元组
        """
        self.logger.info(f"识别Excel列类型")
        
        # 使用正则表达式识别列类型
        cation_pattern = re.compile(r'^cation\d+$', re.IGNORECASE)
        anion_pattern = re.compile(r'^anion\d+$', re.IGNORECASE)
        sol_pattern = re.compile(r'^sol\d+$', re.IGNORECASE)
        
        cation_cols = [header for header in headers if cation_pattern.match(header)]
        anion_cols = [header for header in headers if anion_pattern.match(header)]
        sol_cols = [header for header in headers if sol_pattern.match(header)]
        name_col = 'name' if 'name' in headers else None
        
        self.logger.debug(f"识别到的阳离子列: {cation_cols}")
        self.logger.debug(f"识别到的阴离子列: {anion_cols}")
        self.logger.debug(f"识别到的溶剂列: {sol_cols}")
        
        if not cation_cols:
            self.logger.warning("Excel表中未找到阳离子列")
        if not anion_cols:
            self.logger.warning("Excel表中未找到阴离子列")
        if not name_col:
            self.logger.warning("Excel表中未找到name列")
        
        return cation_cols, anion_cols, sol_cols, name_col
    
    def extract_simulation_parameters(self, row: pd.Series) -> Dict[str, Any]:
        """
        从Excel行中提取模拟参数
        
        Args:
            row: Excel数据行
            
        Returns:
            包含模拟参数的字典
        """
        params = {}
        
        # 提取常见的模拟参数
        param_keys = [
            'box_size', 'box_size_x', 'box_size_y', 'box_size_z',  # 盒子尺寸
            'concentration', 'temperature',                         # 浓度和温度
            'solvent_ratio', 'cation_anion_ratio',                  # 比例参数
            'density', 'cutoff'                                     # 其他参数
        ]
        
        for key in param_keys:
            if key in row and not pd.isna(row[key]):
                params[key] = row[key]
        
        # 处理盒子尺寸 - 支持单一box_size或三维尺寸
        if 'box_size' in params:
            # 单一盒子尺寸（立方体）
            box_size = params['box_size']
            params['box_size'] = (box_size, box_size, box_size)
        elif all(f'box_size_{dim}' in params for dim in ['x', 'y', 'z']):
            # 三维盒子尺寸
            params['box_size'] = (
                params['box_size_x'],
                params['box_size_y'],
                params['box_size_z']
            )
        
        # 提取组分比例参数
        ratio_params = {}
        for col in row.index:
            # 匹配像 cation1_ratio, anion2_ratio, sol3_ratio 这样的列
            ratio_match = re.match(r'(cation\d+|anion\d+|sol\d+)_ratio', col)
            if ratio_match and not pd.isna(row[col]):
                component_name = ratio_match.group(1)
                ratio_params[component_name] = row[col]
        
        if ratio_params:
            params['component_ratios'] = ratio_params
        
        self.logger.debug(f"提取的模拟参数: {params}")
        return params
    
    def calculate_molecules(self, row: pd.Series) -> Dict[str, Dict[str, int]]:
        """
        根据Excel行数据计算所需的分子数量
        
        Args:
            row: Excel数据行
            
        Returns:
            包含计算出的分子数量的字典
        """
        self.logger.info("计算分子数量")
        
        # 提取模拟参数
        params = self.extract_simulation_parameters(row)
        
        # 检查是否已经提供了分子数量
        has_existing_counts = False
        cation_cols = [col for col in row.index if col.lower().startswith('cation') and not col.endswith('_ratio') 
                      and not col.endswith('_name') and not col.endswith('_smile')]
        anion_cols = [col for col in row.index if col.lower().startswith('anion') and not col.endswith('_ratio')
                     and not col.endswith('_name') and not col.endswith('_smile')]
        sol_cols = [col for col in row.index if col.lower().startswith('sol') and not col.endswith('_ratio')
                   and not col.endswith('_name') and not col.endswith('_smile')]
        
        for col in cation_cols + anion_cols + sol_cols:
            if isinstance(row[col], (int, float)) and not pd.isna(row[col]):
                has_existing_counts = True
                break
        
        if has_existing_counts:
            self.logger.info("使用Excel中提供的分子数量")
            molecule_counts = {
                'cations': {col: int(row[col]) for col in cation_cols if isinstance(row[col], (int, float)) and not pd.isna(row[col])},
                'anions': {col: int(row[col]) for col in anion_cols if isinstance(row[col], (int, float)) and not pd.isna(row[col])},
                'solvents': {col: int(row[col]) for col in sol_cols if isinstance(row[col], (int, float)) and not pd.isna(row[col])}
            }
            return molecule_counts
        
        # 检测是否存在比例计算模式（基于主阳离子）
        if 'component_ratios' in params and 'concentration' in params and 'box_size' in params:
            # 使用基于主阳离子的比例计算
            self.logger.info("使用基于主阳离子的比例计算")
            
            # 检查是否有单一盒子尺寸
            box_size = params['box_size']
            if isinstance(box_size, tuple) and all(x == box_size[0] for x in box_size):
                # 三个维度相等，是立方体
                box_size = box_size[0]
            else:
                # 不是立方体，取平均值作为立方体边长
                box_size = (box_size[0] * box_size[1] * box_size[2]) ** (1/3)
                self.logger.warning(f"盒子不是立方体，使用等效立方体边长: {box_size:.2f}")
                
            # 默认cation1的比例为1.0
            component_ratios = params['component_ratios']
            if 'cation1' not in component_ratios:
                component_ratios['cation1'] = 1.0
            
            # 调用分子计算器
            result = self.molecule_calculator.calculate_from_primary_cation(
                box_size=box_size,
                primary_cation_concentration=params['concentration'],
                component_ratios=component_ratios,
                length_unit='angstrom',
                concentration_unit='mol/L'
            )
            
            self.logger.info(f"基于主阳离子的计算结果: {result}")
            return result
        
        # 根据模拟参数自动计算分子数量
        elif 'box_size' in params and 'concentration' in params:
            # 准备计算所需的参数
            box_size = params['box_size']
            concentration = params['concentration']
            solvent_ratio = params.get('solvent_ratio', 0.0)
            cation_anion_ratio = params.get('cation_anion_ratio', 1.0)
            
            # 调用分子计算器
            if solvent_ratio == 0.0 and sol_cols:
                # 如果没有提供溶剂比例但有溶剂列，默认使用一个标准比例
                solvent_ratio = 10.0
                self.logger.warning(f"未提供溶剂比例，使用默认值: {solvent_ratio}")
            
            # 创建溶剂数据
            solvent_data = {}
            for sol_col in sol_cols:
                # 检查是否有专门的溶剂比例
                sol_ratio_col = f"{sol_col}_ratio"
                if sol_ratio_col in row and not pd.isna(row[sol_ratio_col]):
                    # 使用指定的溶剂比例
                    solvent_data[sol_col] = row[sol_ratio_col]
                else:
                    # 如果Excel中没有提供具体比例，假设所有溶剂比例相等
                    solvent_data[sol_col] = 1.0 / len(sol_cols) if len(sol_cols) > 0 else 0.0
            
            # 计算分子数量
            result = self.molecule_calculator.calculate_from_system_definition(
                box_size=box_size,
                cation_concentration=concentration,
                solvent_data=solvent_data,
                cation_types=len(cation_cols),
                anion_types=len(anion_cols),
                length_unit='angstrom',
                concentration_unit='mol/L'
            )
            
            # 映射计算结果到列名
            molecule_counts = {
                'cations': {},
                'anions': {},
                'solvents': {}
            }
            
            # 映射阳离子
            for i, col in enumerate(cation_cols):
                key = f'cation{i+1}'
                if key in result['cations']:
                    molecule_counts['cations'][col] = result['cations'][key]
            
            # 映射阴离子
            for i, col in enumerate(anion_cols):
                key = f'anion{i+1}'
                if key in result['anions']:
                    molecule_counts['anions'][col] = result['anions'][key]
            
            # 映射溶剂
            if solvent_data:
                for i, col in enumerate(sol_cols):
                    solvent_key = list(result['solvents'].keys())[i] if i < len(result['solvents']) else None
                    if solvent_key:
                        molecule_counts['solvents'][col] = result['solvents'][solvent_key]
            
            self.logger.info(f"计算得到的分子数量: {molecule_counts}")
            return molecule_counts
        else:
            self.logger.error("无法计算分子数量，缺少必要参数")
            return {
                'cations': {},
                'anions': {},
                'solvents': {}
            }
    
    def generate_new_row_with_counts(self, row: pd.Series) -> pd.Series:
        """
        生成包含计算出的分子数量的新行
        
        Args:
            row: 原始Excel数据行
            
        Returns:
            包含计算出的分子数量的新行
        """
        self.logger.info("生成带有计算数量的新行")
        
        # 创建新行，复制原始行的所有数据
        new_row = row.copy()
        
        # 计算分子数量
        molecule_counts = self.calculate_molecules(row)
        
        # 更新新行中的分子数量
        for category, counts in molecule_counts.items():
            for column, count in counts.items():
                new_row[column] = count
        
        self.logger.debug(f"生成的新行: {new_row}")
        return new_row
    
    def process_excel_file(self, file_path: str, sheet_name: Optional[str] = None, 
                         header_row: int = 1, output_file: Optional[str] = None) -> pd.DataFrame:
        """
        处理Excel文件，自动计算分子数量
        
        Args:
            file_path: Excel文件路径
            sheet_name: 工作表名称，如果为None则使用第一个工作表
            header_row: 表头行号（从1开始）
            output_file: 输出文件路径，如果提供则保存处理后的Excel
            
        Returns:
            处理后的DataFrame
        """
        self.logger.info(f"处理Excel文件: {file_path}")
        
        # 读取Excel文件
        df = self.read_excel(file_path, sheet_name, header_row)
        
        # 处理每一行，计算分子数量
        new_rows = []
        for _, row in df.iterrows():
            new_row = self.generate_new_row_with_counts(row)
            new_rows.append(new_row)
        
        # 创建新的DataFrame
        new_df = pd.DataFrame(new_rows)
        
        # 如果提供了输出文件路径，保存结果
        if output_file:
            new_df.to_excel(output_file, index=False)
            self.logger.info(f"处理后的Excel已保存到: {output_file}")
        
        return new_df 