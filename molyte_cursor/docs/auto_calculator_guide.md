# Molyte自动分子计算器使用指南

本指南介绍如何使用Molyte的自动分子计算功能，该功能可以根据物理化学参数（如体积、浓度、比例等）自动计算模拟所需的分子数量，无需在Excel中手动计算。

## 目录

- [功能简介](#功能简介)
- [基本概念](#基本概念)
- [Excel模板格式](#excel模板格式)
  - [必填字段](#必填字段)
  - [组分信息字段](#组分信息字段)
  - [示例模板](#示例模板)
- [使用方法](#使用方法)
  - [使用ExcelReader](#使用excelreader)
  - [使用ExcelProcessor](#使用excelprocessor)
  - [直接使用MoleculeCalculator](#直接使用moleculecalculator)
- [示例](#示例)
- [常见问题](#常见问题)

## 功能简介

传统上，用户需要在Excel文件中手动计算分子动力学模拟所需的各种分子数量。这不仅繁琐，还容易出错。Molyte的自动分子计算功能可以：

1. 根据指定的体积、浓度和比例自动计算分子数量
2. 支持多种阳离子、阴离子和溶剂的配比计算
3. 基于主阳离子和比例计算各组分数量
4. 自动处理Excel输入，无需用户手动填写分子数量
5. 保留用户手动指定数量的灵活性

## 基本概念

### 计算原理

分子数量的计算基于以下基本原理：

- **主阳离子数量** = 浓度 × 体积 × 阿伏伽德罗常数
- **其他组分数量** = 主阳离子数量 × 相应的比例

其中：
- 浓度单位为mol/L（或mol/m³）
- 体积单位为立方埃（立方体盒子）
- 阿伏伽德罗常数为6.02214076×10²³ mol⁻¹

## Excel模板格式

### 必填字段

Excel模板中必须提供以下字段：

1. **name**：系统名称（如"System1"、"Li-TFSI-EC-DMC"等）
2. **temperature**：模拟温度（单位：K）
3. **box_size**：立方体盒子边长（单位：埃，默认值：40，有效范围：0-100）
4. **concentration**：第一种阳离子(cation1)的浓度（单位：mol/L）

### 组分信息字段

对于每种组分（阳离子、阴离子、溶剂），需要提供以下信息：

1. **[组分]_name**：组分名称（例如：cation1_name="Li"）
2. **[组分]_smile**：组分的SMILE字符串（例如：cation1_smile="[Li+]"）
3. **[组分]_ratio**：相对于第一种阳离子的比例（例如：cation1_ratio=1.0）

**注意**：第一种阳离子(cation1)的比例必须设为1.0，其他所有组分的比例都是相对于第一种阳离子的。

### 示例模板

一个完整的Excel模板示例：

| name | temperature | box_size | concentration | cation1_name | cation1_smile | cation1_ratio | cation2_name | cation2_smile | cation2_ratio | anion1_name | anion1_smile | anion1_ratio | sol1_name | sol1_smile | sol1_ratio | sol2_name | sol2_smile | sol2_ratio |
|------|-------------|----------|---------------|--------------|---------------|---------------|--------------|---------------|---------------|-------------|-------------|--------------|-----------|------------|------------|-----------|------------|------------|
| Li-TFSI-EC-DMC | 300 | 50 | 1.0 | Li | [Li+] | 1.0 | Na | [Na+] | 0.5 | TFSI | N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F | 1.5 | EC | C1COC(=O)O1 | 7.0 | DMC | COC(=O)OC | 3.0 |

**计算过程**：
1. 系统先验证box_size是否在有效范围内（0-100埃），若不在范围内则拒绝计算
2. 系统先计算第一种阳离子(Li)的数量，基于浓度(1.0 mol/L)和盒子体积(50×50×50 Å³)
3. 然后根据比例计算其他组分的数量：
   - Na数量 = Li数量 × 0.5
   - TFSI数量 = Li数量 × 1.5
   - EC数量 = Li数量 × 7.0
   - DMC数量 = Li数量 × 3.0

## 使用方法

### 使用ExcelReader

`ExcelReader`类已经集成了自动计算功能，可以直接使用：

```python
from molyte_cursor.src.io.excel_reader import ExcelReader

# 初始化Excel读取器，启用自动计算
reader = ExcelReader(auto_calculate=True)

# 读取Excel文件
systems = reader.read_excel('input.xlsx')

# 使用读取的系统
for system in systems:
    print(f"系统: {system.name}")
    print(f"阳离子: {[f'{c.name}({c.count})' for c in system.cations]}")
    print(f"阴离子: {[f'{a.name}({a.count})' for a in system.anions]}")
    print(f"溶剂: {[f'{s.name}({s.count})' for s in system.solvents]}")
```

### 使用ExcelProcessor

如果您想预处理Excel文件，将自动计算的结果保存到新的Excel文件中：

```python
from molyte_cursor.src.io.excel_processor import ExcelProcessor

# 初始化Excel处理器
processor = ExcelProcessor()

# 处理Excel文件，生成带有计算结果的新Excel
result_df = processor.process_excel_file(
    file_path='input.xlsx',
    output_file='output.xlsx'
)
```

### 直接使用MoleculeCalculator

对于更复杂的场景，可以直接使用`MoleculeCalculator`类：

```python
from molyte_cursor.src.utils.molecule_calculator import MoleculeCalculator

# 初始化分子计算器
calculator = MoleculeCalculator()

# 基于主阳离子和比例的计算
result = calculator.calculate_from_primary_cation(
    box_size=50,                      # 立方体盒子边长（埃），必须在0-100范围内
    primary_cation_concentration=1.0, # 主阳离子浓度（mol/L）
    component_ratios={                # 各组分相对于主阳离子的比例
        'cation1': 1.0,               # 主阳离子(Li)，比例为1.0
        'cation2': 0.5,               # 第二种阳离子(Na)
        'anion1': 1.5,                # 第一种阴离子(TFSI)
        'sol1': 7.0,                  # 第一种溶剂(EC)
        'sol2': 3.0                   # 第二种溶剂(DMC)
    }
)

print(f"阳离子: {result['cations']}")
print(f"阴离子: {result['anions']}")
print(f"溶剂: {result['solvents']}")
```

## 示例

完整示例可以参考以下文件：

```bash
cd molyte_cursor
python examples/demo_ratio_based_calculator.py
```

该脚本演示了：
1. 直接使用`MoleculeCalculator`进行计算
2. 创建和处理演示Excel文件
3. 使用`ExcelProcessor`处理Excel文件
4. 使用`ExcelReader`读取Excel并应用自动计算

## 常见问题

### 1. 自动计算与手动数值

**问题**：如果Excel中已有分子数量，自动计算会覆盖这些值吗？

**回答**：不会。系统会优先使用Excel中已有的数值。只有当相关列为空时，才会进行自动计算。

### 2. 支持的单位

**问题**：自动计算支持哪些单位？

**回答**：
- 体积：立方埃(angstrom3)
- 浓度：mol/L(摩尔/升)
- 温度：K(开尔文)

### 3. 多种溶剂的比例设置

**问题**：如何设置多种溶剂的比例？

**回答**：对每种溶剂设置`sol[n]_ratio`字段，表示该溶剂与第一种阳离子的比例。例如，`sol1_ratio=7.0`表示第一种溶剂与第一种阳离子的比例为7:1。

### 4. SMILE字符串格式

**问题**：SMILE字符串需要遵循什么格式？

**回答**：需要使用标准的SMILE格式。例如，Li+表示为"[Li+]"，TFSI表示为"N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F"。

### 5. 盒子尺寸限制

**问题**：为什么有盒子尺寸的限制？

**回答**：盒子尺寸(box_size)有以下限制：
- 默认值：如果不指定，系统会使用40埃作为默认值
- 有效范围：盒子尺寸必须大于0且小于或等于100埃
- 限制原因：
  1. 过小的盒子可能无法容纳所有指定的分子
  2. 过大的盒子会导致计算资源消耗过多，并可能生成不合理的分子数量
  3. 程序设计上限制在100埃以内，以确保计算结果的合理性和计算效率

如果遇到盒子尺寸限制的错误，请调整box_size参数至合理范围。 