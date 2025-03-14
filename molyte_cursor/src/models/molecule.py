"""
分子模型模块，定义分子相关的数据结构
"""
from dataclasses import dataclass
from typing import List, Dict, Optional, Any
from pathlib import Path

@dataclass
class Molecule:
    """分子基本类"""
    name: str                  # 分子名称
    original_name: str         # 原始名称（可能与name相同）
    number: int                # 分子数量
    ratio: Optional[float] = None  # 比例（如果有的话）
    smile: Optional[str] = None    # SMILE表示（如果有的话）
    
    @classmethod
    def from_dict(cls, data):
        """从字典创建分子实例
        
        Args:
            data: 分子数据字典
            
        Returns:
            分子实例
        """
        return cls(
            name=data.get('name'),
            original_name=data.get('original_name', data.get('name')),
            number=data.get('number', 0),
            ratio=data.get('ratio'),
            smile=data.get('smile')
        )
    
    def to_dict(self):
        """转换为字典
        
        Returns:
            分子数据字典
        """
        result = {
            'name': self.name,
            'original_name': self.original_name,
            'number': self.number
        }
        
        if self.ratio is not None:
            result['ratio'] = self.ratio
            
        if self.smile is not None:
            result['smile'] = self.smile
            
        return result
    
    def is_valid(self):
        """检查分子是否有效
        
        Returns:
            布尔值，True表示有效
        """
        return self.name and self.number > 0

@dataclass
class Cation(Molecule):
    """阳离子类"""
    type: str = "cation"
    
@dataclass
class Anion(Molecule):
    """阴离子类"""
    type: str = "anion"
    
@dataclass
class Solvent(Molecule):
    """溶剂类"""
    type: str = "solvent"

@dataclass
class MoleculeSystem:
    """分子系统类，表示一个完整的模拟系统"""
    name: str                    # 系统名称
    cations: List[Cation]        # 阳离子列表
    anions: List[Anion]          # 阴离子列表
    solvents: List[Solvent]      # 溶剂列表
    temperature: float           # 温度
    box_size: float              # 盒子大小
    additional_params: Dict[str, Any] = None  # 额外参数
    
    def __post_init__(self):
        """初始化后的处理"""
        if self.additional_params is None:
            self.additional_params = {}
    
    def has_valid_cation(self):
        """检查是否有有效的阳离子
        
        Returns:
            布尔值，True表示有有效阳离子
        """
        return any(cation.is_valid() for cation in self.cations)
    
    def to_dict(self):
        """转换为字典
        
        Returns:
            系统数据字典
        """
        return {
            'name': self.name,
            'cations': [cation.to_dict() for cation in self.cations],
            'anions': [anion.to_dict() for anion in self.anions],
            'solvents': [solvent.to_dict() for solvent in self.solvents],
            'temperature': self.temperature,
            'box_size': self.box_size,
            'additional_params': self.additional_params
        }
    
    @classmethod
    def from_dict(cls, data):
        """从字典创建系统实例
        
        Args:
            data: 系统数据字典
            
        Returns:
            系统实例
        """
        return cls(
            name=data.get('name'),
            cations=[Cation.from_dict(c) for c in data.get('cations', [])],
            anions=[Anion.from_dict(a) for a in data.get('anions', [])],
            solvents=[Solvent.from_dict(s) for s in data.get('solvents', [])],
            temperature=data.get('temperature'),
            box_size=data.get('box_size'),
            additional_params=data.get('additional_params', {})
        ) 