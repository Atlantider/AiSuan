# Molyte分子性质计算指南

本文档介绍如何使用Molyte中的高斯处理工具计算分子电子性质（如HOMO、LUMO、ESP等）并进行可视化。

## 目录

- [功能概述](#功能概述)
- [准备工作](#准备工作)
- [命令行工具](#命令行工具)
- [编程接口](#编程接口)
- [结果分析](#结果分析)
- [常见问题](#常见问题)

## 功能概述

Molyte的高斯处理工具提供以下功能：

1. **自动生成输入文件**：从Excel表格中读取溶剂分子信息，自动生成高斯计算输入文件
2. **批量处理**：为多个分子生成批处理脚本，便于批量计算
3. **解析结果**：自动提取高斯输出文件中的关键信息，如HOMO、LUMO能级、能隙、偶极矩等
4. **可视化结果**：生成直观的图表，帮助分析分子性质
5. **系统集成**：与Molyte其他模块紧密集成，实现从分子定义到性质计算的无缝流程

## 准备工作

### 前提条件

使用此功能需要以下条件：

1. **高斯软件**：计算需要安装Gaussian 16或09，用于量子化学计算
2. **Python依赖**：需要安装以下Python库：
   - numpy
   - matplotlib
   - pandas
   - rdkit (用于从SMILE生成3D结构)

### 输入准备

程序需要一个Excel文件，包含溶剂分子信息。Excel文件格式应符合Molyte标准，特别是：

- 每个溶剂必须有`name`和`smile`字段
- 分子结构应通过SMILE字符串表示

## 命令行工具

### 基本用法

```bash
# 从Excel生成高斯输入文件
python -m molyte_cursor.src.analysis.calc_properties path/to/excel_file.xlsx

# 生成输入文件并立即运行高斯计算
python -m molyte_cursor.src.analysis.calc_properties path/to/excel_file.xlsx --run

# 只解析已有的高斯输出文件并可视化
python -m molyte_cursor.src.analysis.calc_properties path/to/excel_file.xlsx --parse-only
```

### 命令行选项

```
-o, --output-dir    指定输出目录，默认为"gaussian_calc"
-r, --run           生成输入文件后立即运行高斯计算
-g, --gaussian      高斯可执行文件路径或命令，默认为"g16"
-m, --method        计算方法，默认为"B3LYP"
-b, --basis         基组，默认为"6-31G(d)"
-c, --calc-type     计算类型，默认为"opt freq"
-v, --verbose       输出详细日志信息
-p, --parse-only    只解析已有的高斯输出文件而不生成新的输入
```

### 使用示例

1. 生成标准计算的输入文件:
   ```bash
   python -m molyte_cursor.src.analysis.calc_properties mysolvent.xlsx -o solvent_calc
   ```

2. 使用特定的计算参数:
   ```bash
   python -m molyte_cursor.src.analysis.calc_properties mysolvent.xlsx -m "wB97XD" -b "6-311G(d,p)" -c "opt freq pop=ESP"
   ```

3. 运行计算(如果已安装高斯):
   ```bash
   python -m molyte_cursor.src.analysis.calc_properties mysolvent.xlsx --run
   ```

4. 解析之前计算的结果:
   ```bash
   python -m molyte_cursor.src.analysis.calc_properties mysolvent.xlsx --parse-only -o previous_calc
   ```

## 编程接口

如果需要在Python脚本中使用该功能，可以直接导入`GaussianProcessor`类:

```python
from molyte_cursor.src.analysis.gaussian_processor import GaussianProcessor

# 创建高斯处理器
processor = GaussianProcessor()

# 从Excel生成高斯输入文件
input_files = processor.read_excel_and_generate_inputs(
    "solvents.xlsx", 
    "gaussian_calc",
    method="B3LYP",
    basis_set="6-31G(d)"
)

# 生成批处理脚本
bash_script, bat_script = processor.generate_batch_script(
    input_files,
    "gaussian_calc"
)

# 解析高斯输出文件
results, summary_df = processor.parse_all_outputs("gaussian_calc")

# 可视化结果
plot_files = processor.visualize_results(results)

# 或者执行完整的流程
summary_df, plot_files = processor.calculate_and_visualize_properties(
    "solvents.xlsx",
    "gaussian_calc",
    run_gaussian=False
)
```

## 结果分析

### 生成的文件

程序会在指定的输出目录生成以下文件：

1. **高斯输入文件**：`.gjf`文件，每个溶剂一个
2. **批处理脚本**：`run_gaussian.sh`(Linux/Mac)和`run_gaussian.bat`(Windows)
3. **结果摘要**：`gaussian_results_summary.csv`，包含各分子的计算结果
4. **可视化图表**：
   - `homo_lumo_comparison.png`: HOMO-LUMO能级对比图
   - `energy_gap_comparison.png`: 能隙对比图
   - `dipole_moment_comparison.png`: 偶极矩对比图
   - `homo_dipole_correlation.png`: HOMO能级与偶极矩相关性图

### 结果解读

1. **HOMO-LUMO能级**：
   - HOMO (最高占据分子轨道) 能级越高，表示分子越容易给出电子（还原能力越强）
   - LUMO (最低未占据分子轨道) 能级越低，表示分子越容易接受电子（氧化能力越强）

2. **能隙**：
   - 能隙 = LUMO - HOMO，表示分子的稳定性和反应活性
   - 较小的能隙意味着较高的反应活性和较低的化学稳定性

3. **偶极矩**：
   - 表示电荷分布的不均匀性，与分子的极性相关
   - 较大的偶极矩通常意味着更好的溶解性（在极性溶剂中）

## 常见问题

### 1. 高斯计算失败

问题: 高斯计算失败，没有正常终止。

回答: 可能的原因包括:
- 初始结构不合理，导致SCF不收敛
- 内存或磁盘空间不足
- 高斯许可证问题

解决方法:
- 检查高斯错误信息
- 尝试使用不同的初始结构
- 调整计算参数，如增加SCF迭代次数

### 2. 找不到特定分子轨道能级

问题: 解析结果中缺少HOMO、LUMO等数据。

回答: 
- 确保计算类型包含足够的输出信息（如`pop=full`）
- 检查高斯输出文件是否完整
- 对于开壳层体系，程序目前只提取Alpha轨道信息

### 3. RDKit导入错误

问题: 生成输入文件时出现RDKit导入错误。

回答:
- 安装RDKit: `conda install -c conda-forge rdkit`
- 如不需要自动生成3D结构，可以手动准备分子结构文件

### 4. SMILE解析问题

问题: 某些SMILE字符串无法被正确解析。

回答:
- 确保SMILE格式正确，特别是对于复杂或带电的分子
- 可以使用在线工具验证SMILE字符串（如ChemDraw）
- 对于有问题的分子，考虑手动准备结构文件 