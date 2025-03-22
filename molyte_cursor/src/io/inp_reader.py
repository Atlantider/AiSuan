"""
INP文件读取器

用于解析Web前端生成的INP格式文件，转换为molyte-cursor可以处理的配置数据。
"""

import os
import re
import logging
from typing import Dict, List, Any, Optional, Tuple
import traceback

class INPReader:
    """
    读取和解析INP文件的类
    
    INP文件是一种简单的文本格式，用于描述电解液配方和模拟参数。
    格式为：
    ```
    START
    name 配方名称
    T 298
    Box_size 40
    concentration 1.0
    Cation1_name Li
    Cation1_ratio 1.0
    Anion1_name PF6
    Anion1_ratio 1.0
    Sol1_name EC
    Sol1_smile C1COC(=O)O1
    Sol1_ratio 5.0
    ...
    END
    ```
    """
    
    def __init__(self, file_path: str = None, content: str = None):
        """
        初始化INP读取器
        
        Args:
            file_path: INP文件路径，与content参数二选一
            content: INP文件内容字符串，与file_path参数二选一
        """
        self.file_path = file_path
        self.content = content
        self.logger = logging.getLogger(__name__)
        
        if not file_path and not content:
            raise ValueError("必须提供file_path或content参数中的一个")
    
    def read(self) -> Dict[str, Any]:
        """
        读取并解析INP文件内容
        
        Returns:
            Dict: 包含解析后参数的字典
        """
        # 如果提供了content，直接使用；否则从文件读取
        if self.content:
            content = self.content
        else:
            with open(self.file_path, 'r', encoding='utf-8') as f:
                content = f.read()
        
        return self.parse_content(content)
    
    def parse_content(self, content: str) -> Dict[str, Any]:
        """
        解析INP文件内容字符串
        
        Args:
            content: INP文件内容
            
        Returns:
            Dict: 包含解析后参数的字典
        """
        # 初始化结果字典
        result = {
            'general': {},
            'cations': [],
            'anions': [],
            'solvents': []
        }
        
        # 检查内容格式
        if 'START' not in content or 'END' not in content:
            self.logger.warning("INP文件缺少START或END标记")
        
        # 提取START和END之间的内容
        pattern = r'START\s*(.*?)\s*END'
        match = re.search(pattern, content, re.DOTALL)
        if match:
            content = match.group(1)
        
        # 解析每一行
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            try:
                # 先检查是否是特殊的SMILE行
                smile_match = re.search(r'^(sol\d+_smile)\s*(?:=\s*)?(.+)$', line, re.IGNORECASE)
                if smile_match:
                    key = smile_match.group(1).lower()
                    value = smile_match.group(2).strip()
                    
                    solvent_index = int(re.search(r'sol(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['solvents']) <= solvent_index:
                        result['solvents'].append({})
                    
                    # 检查对应的溶剂名称
                    solvent_name = ""
                    if solvent_index < len(result['solvents']):
                        solvent_name = result['solvents'][solvent_index].get('name', '').upper()
                    
                    result['solvents'][solvent_index]['smile'] = value
                    
                    self.logger.debug(f"解析SMILE: {key}={value}, 溶剂名称={solvent_name}")
                    continue
                
                # 处理其他类型的行
                if '=' in line:
                    # 处理"key = value"格式
                    parts = line.split('=', 1)
                    key = parts[0].strip().lower()
                    value = parts[1].strip()
                else:
                    # 处理原始的"key value"格式
                    parts = line.split(' ', 1)
                    key = parts[0].strip().lower()
                    value = parts[1].strip() if len(parts) > 1 else ""
                
                # 解析阳离子
                if re.search(r'cation\d+_name', key, re.IGNORECASE):
                    cation_index = int(re.search(r'cation(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['cations']) <= cation_index:
                        result['cations'].append({})
                    result['cations'][cation_index]['name'] = value
                
                elif re.search(r'cation\d+_ratio', key, re.IGNORECASE):
                    cation_index = int(re.search(r'cation(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['cations']) <= cation_index:
                        result['cations'].append({})
                    result['cations'][cation_index]['ratio'] = float(value)
                
                # 解析阴离子
                elif re.search(r'anion\d+_name', key, re.IGNORECASE):
                    anion_index = int(re.search(r'anion(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['anions']) <= anion_index:
                        result['anions'].append({})
                    result['anions'][anion_index]['name'] = value
                
                elif re.search(r'anion\d+_ratio', key, re.IGNORECASE):
                    anion_index = int(re.search(r'anion(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['anions']) <= anion_index:
                        result['anions'].append({})
                    result['anions'][anion_index]['ratio'] = float(value)
                
                # 解析溶剂名称
                elif re.search(r'sol\d+_name', key, re.IGNORECASE):
                    solvent_index = int(re.search(r'sol(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['solvents']) <= solvent_index:
                        result['solvents'].append({})
                    result['solvents'][solvent_index]['name'] = value
                
                # 解析溶剂比例
                elif re.search(r'sol\d+_ratio', key, re.IGNORECASE):
                    solvent_index = int(re.search(r'sol(\d+)', key, re.IGNORECASE).group(1)) - 1
                    while len(result['solvents']) <= solvent_index:
                        result['solvents'].append({})
                    result['solvents'][solvent_index]['ratio'] = float(value)
                
                # 常规参数
                else:
                    # 尝试转换数值类型
                    try:
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except ValueError:
                        pass
                    
                    result['general'][key] = value
            
            except Exception as e:
                self.logger.warning(f"解析行'{line}'时出错: {str(e)}")
                self.logger.debug(f"异常详情: {traceback.format_exc()}")
        
        # 计算浓度（如果没有直接提供）
        if 'concentration' not in result['general'] and len(result['cations']) > 0:
            # 假设第一个阳离子的浓度为基准
            result['general']['concentration'] = 1.0
        
        return result
    
    def convert_to_molyte_config(self) -> Dict[str, Any]:
        """
        将INP解析结果转换为molyte-cursor所需的配置格式
        
        Returns:
            Dict: molyte-cursor配置字典
        """
        # 读取并解析INP文件
        inp_data = self.read()
        
        # 构建配置
        general = inp_data['general']
        
        # 基础配置
        config = {
            'formulation_name': general.get('name', 'Unnamed Formulation'),
            'temperature': general.get('t', general.get('temperature', 298)),
            'box_size': general.get('box_size', 40),
            'concentration': general.get('concentration', 1.0),
            'time_step': general.get('time_step', 1.0),
            'equilibration_steps': general.get('equilibration_steps', 1000000),
            'production_steps': general.get('production_steps', 2000000),
            'cutoff': general.get('cutoff', 12.0),
            'pressure': general.get('pressure', 1.0),
            'salts': [],
            'solvents': []
        }
        
        # 确保基本参数都有值
        base_concentration = config['concentration']
        
        # 处理阳离子和阴离子，组合为盐
        for i, cation in enumerate(inp_data['cations']):
            if not cation.get('name'):
                continue
                
            # 找到对应的阴离子
            anion = inp_data['anions'][i] if i < len(inp_data['anions']) else {'name': 'PF6', 'ratio': 1.0}
            
            # 计算浓度
            cation_ratio = cation.get('ratio', 1.0)
            salt_concentration = base_concentration * cation_ratio
            
            # 移除离子名称中的+和-符号
            cation_name = cation['name'].replace('+', '').replace('-', '')
            anion_name = anion['name'].replace('+', '').replace('-', '')
            
            # 创建盐配置
            salt = {
                'name': f"{cation_name}{anion_name}",
                'cation': cation_name,
                'anion': anion_name,
                'concentration': salt_concentration
            }
            
            config['salts'].append(salt)
            
            # 展平cation和anion信息到主配置中 (兼容ElectrolyteWorkflow的访问方式)
            idx = i + 1  # 1-indexed
            config[f"Cation{idx}_name"] = cation_name
            config[f"Cation{idx}_ratio"] = cation.get('ratio', 1.0)
            if 'smile' in cation:
                config[f"Cation{idx}_smile"] = cation['smile']
                
            config[f"Anion{idx}_name"] = anion_name
            config[f"Anion{idx}_ratio"] = anion.get('ratio', 1.0)
            if 'smile' in anion:
                config[f"Anion{idx}_smile"] = anion['smile']
        
        # 处理溶剂
        for i, solvent in enumerate(inp_data['solvents']):
            if not solvent.get('name'):
                continue
                
            # 获取溶剂浓度比例
            ratio = solvent.get('ratio', 5.0)
            
            # 创建溶剂配置
            solvent_config = {
                'name': solvent['name'],
                'smile': solvent.get('smile', ''),
                'concentration': ratio
            }
            
            config['solvents'].append(solvent_config)
            
            # 展平solvent信息到主配置中 (兼容ElectrolyteWorkflow的访问方式)
            idx = i + 1  # 1-indexed
            config[f"Sol{idx}_name"] = solvent['name']
            config[f"Sol{idx}_ratio"] = solvent.get('ratio', 5.0)
            if 'smile' in solvent:
                config[f"Sol{idx}_smile"] = solvent['smile']
        
        # 如果没有溶剂，添加默认溶剂
        if not config['solvents']:
            config['solvents'].append({
                'name': 'EC',
                'smile': 'C1COC(=O)O1',
                'concentration': 5.0
            })
            # 同时展平到主配置
            config["Sol1_name"] = 'EC'
            config["Sol1_smile"] = 'C1COC(=O)O1'
            config["Sol1_ratio"] = 5.0
        else:
            # 验证所有溶剂的SMILE字符串
            for i, solvent in enumerate(config['solvents']):
                # 确保SMILE字符串存在
                if 'smile' not in solvent or not solvent['smile']:
                    self.logger.warning(f"溶剂 {solvent.get('name', f'未知溶剂{i}')} 缺少SMILE字符串，将尝试提供默认值")
                    
                    # 检查是否是已知溶剂
                    if solvent.get('name', '').upper() == 'EC':
                        solvent['smile'] = 'C1COC(=O)O1'
                        self.logger.info(f"为EC溶剂添加默认SMILE字符串: 'C1COC(=O)O1'")
                    elif solvent.get('name', '').upper() == 'DMC':
                        solvent['smile'] = 'COC(=O)OC'
                        self.logger.info(f"为DMC溶剂添加默认SMILE字符串: 'COC(=O)OC'")
                    elif solvent.get('name', '').upper() == 'DEC':
                        solvent['smile'] = 'CCOC(=O)OCC'
                        self.logger.info(f"为DEC溶剂添加默认SMILE字符串: 'CCOC(=O)OCC'")
                    elif solvent.get('name', '').upper() == 'EMC':
                        solvent['smile'] = 'CCOC(=O)OC'
                        self.logger.info(f"为EMC溶剂添加默认SMILE字符串: 'CCOC(=O)OC'")
                    elif solvent.get('name', '').upper() == 'PC':
                        solvent['smile'] = 'CC1COC(=O)O1'
                        self.logger.info(f"为PC溶剂添加默认SMILE字符串: 'CC1COC(=O)O1'")
                
                # 验证SMILE字符串格式
                if 'smile' in solvent and solvent['smile']:
                    # 检查SMILE字符串中的括号是否配对
                    parenthesis = 0
                    brackets = 0
                    for c in solvent['smile']:
                        if c == '(':
                            parenthesis += 1
                        elif c == ')':
                            parenthesis -= 1
                        elif c == '[':
                            brackets += 1
                        elif c == ']':
                            brackets -= 1
                    
                    if parenthesis != 0 or brackets != 0:
                        self.logger.warning(f"溶剂 {solvent.get('name', f'未知溶剂{i}')} 的SMILE字符串括号不匹配: '{solvent['smile']}'")
                
                # 更新配置字典
                idx = i + 1  # 1-indexed
                config[f"Sol{idx}_name"] = solvent['name']
                config[f"Sol{idx}_ratio"] = solvent.get('ratio', 5.0)
                if 'smile' in solvent:
                    config[f"Sol{idx}_smile"] = solvent['smile']
        
        # 处理计算类型
        if 'calculation_types' in general:
            calculation_types = general['calculation_types']
            if isinstance(calculation_types, str):
                calculation_types = [t.strip() for t in calculation_types.split(',')]
            config['calculation_types'] = calculation_types
        else:
            config['calculation_types'] = ['conductivity']
        
        return config

def read_inp_file(file_path: str) -> Dict[str, Any]:
    """
    读取INP文件并返回molyte-cursor配置
    
    Args:
        file_path: INP文件路径
        
    Returns:
        Dict: molyte-cursor配置字典
    """
    reader = INPReader(file_path=file_path)
    return reader.convert_to_molyte_config()

def parse_inp_content(content: str) -> Dict[str, Any]:
    """
    解析INP内容字符串并返回molyte-cursor配置
    
    Args:
        content: INP文件内容
        
    Returns:
        Dict: molyte-cursor配置字典
    """
    reader = INPReader(content=content)
    return reader.convert_to_molyte_config() 