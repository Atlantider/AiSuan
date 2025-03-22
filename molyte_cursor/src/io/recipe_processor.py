"""
配方处理器模块，用于处理从网页提交的JSON格式电解液配方
"""
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
from ..utils.logger import Logger

@dataclass
class MoleculeComponent:
    """分子组分数据类"""
    name: str
    smile: Optional[str]  # 对于离子可能没有SMILE
    number: int
    ratio: float
    charge: int = 0

@dataclass
class ElectrolyteRecipe:
    """电解液配方数据类"""
    name: str
    temperature: float
    box_size: float
    cations: List[MoleculeComponent]
    anions: List[MoleculeComponent]
    solvents: List[MoleculeComponent]

class RecipeProcessor:
    """配方处理器类，用于解析和验证JSON格式的电解液配方"""
    
    def __init__(self):
        """初始化配方处理器"""
        self.logger = Logger().get_logger()
    
    def load_recipe_from_json(self, json_data: Union[str, Dict]) -> ElectrolyteRecipe:
        """从JSON数据加载电解液配方
        
        Args:
            json_data: JSON字符串或字典格式的配方数据
            
        Returns:
            ElectrolyteRecipe对象
            
        Raises:
            ValueError: 当配方数据格式不正确时
        """
        if isinstance(json_data, str):
            try:
                data = json.loads(json_data)
            except json.JSONDecodeError as e:
                self.logger.error(f"JSON解析错误: {e}")
                raise ValueError(f"无效的JSON格式: {e}")
        else:
            data = json_data
            
        # 验证必需字段
        required_fields = ['name', 'temperature', 'box_size']
        for field in required_fields:
            if field not in data:
                raise ValueError(f"缺少必需字段: {field}")
        
        # 解析分子组分
        cations = self._parse_components(data.get('cations', []))
        anions = self._parse_components(data.get('anions', []))
        solvents = self._parse_components(data.get('solvents', []))
        
        # 验证至少有一个阳离子
        if not cations:
            raise ValueError("配方必须包含至少一个阳离子")
            
        return ElectrolyteRecipe(
            name=data['name'],
            temperature=float(data['temperature']),
            box_size=float(data['box_size']),
            cations=cations,
            anions=anions,
            solvents=solvents
        )
    
    def _parse_components(self, components_data: List[Dict]) -> List[MoleculeComponent]:
        """解析分子组分数据
        
        Args:
            components_data: 组分数据列表
            
        Returns:
            MoleculeComponent对象列表
        """
        components = []
        for comp in components_data:
            try:
                component = MoleculeComponent(
                    name=comp['name'],
                    smile=comp.get('smile'),  # 离子可能没有SMILE
                    number=int(comp['number']),
                    ratio=float(comp.get('ratio', 1.0)),
                    charge=int(comp.get('charge', 0))
                )
                if component.number > 0:
                    components.append(component)
                else:
                    self.logger.warning(f"跳过数量为0的组分: {component.name}")
            except (KeyError, ValueError) as e:
                self.logger.warning(f"解析组分数据时出错: {e}")
                continue
        return components
    
    def validate_recipe(self, recipe: ElectrolyteRecipe) -> bool:
        """验证配方的有效性
        
        Args:
            recipe: 电解液配方对象
            
        Returns:
            bool: 配方是否有效
        """
        # 验证温度范围
        if not (200 <= recipe.temperature <= 500):
            self.logger.error(f"温度 {recipe.temperature}K 超出有效范围 [200K, 500K]")
            return False
            
        # 验证盒子大小
        if not (20 <= recipe.box_size <= 100):
            self.logger.error(f"盒子大小 {recipe.box_size}Å 超出有效范围 [20Å, 100Å]")
            return False
            
        # 验证电荷平衡
        total_charge = sum(cat.number * cat.charge for cat in recipe.cations) + \
                      sum(an.number * an.charge for an in recipe.anions)
        if total_charge != 0:
            self.logger.error(f"系统总电荷不为零: {total_charge}")
            return False
            
        return True 