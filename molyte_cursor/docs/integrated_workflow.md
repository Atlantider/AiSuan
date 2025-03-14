# Molyte集成工作流使用指南

本文档介绍如何使用Molyte集成工作流，该功能将LAMMPS模拟和分子性质计算集成到一个统一的流程中，实现从输入Excel到结果分析的一站式解决方案。

## 目录

- [概述](#概述)
- [工作目录结构](#工作目录结构)
- [预处理阶段](#预处理阶段)
- [后处理阶段](#后处理阶段)
- [使用示例](#使用示例)
- [常见问题](#常见问题)

## 概述

Molyte集成工作流将之前的多个独立功能集成到一个统一的界面中，主要提供两个阶段的工作流：

1. **预处理阶段**：
   - 从Excel读取系统信息
   - 设置并运行LAMMPS模拟
   - 生成高斯输入文件并可选择运行量子化学计算

2. **后处理阶段**：
   - 处理LAMMPS模拟结果，如修复周期性边界问题
   - 分析高斯计算结果，提取电子性质并生成可视化图表
   - 将结果整合到标准的结果目录

集成工作流使用统一的目录结构，确保所有结果文件都被适当组织，便于后续查找和分析。

## 工作目录结构

集成工作流使用以下标准目录结构：

```
工作空间根目录/
├── Lammps_Workspace/  # LAMMPS计算的工作目录
├── Molecule_Workspace/ # 分子计算的工作目录
├── Lammps_Results/    # LAMMPS结果和分析
└── Gaussian_Results/  # 高斯计算结果和分析
```

这些目录会在首次运行时自动创建。您可以通过`--workspace`参数指定工作空间根目录，默认为当前目录。

## 预处理阶段

预处理阶段用于从Excel文件中读取系统信息，设置并可选运行LAMMPS模拟和高斯计算。

### 命令行参数

```
--preprocess          执行预处理阶段
-i, --input           输入Excel文件路径（必需）
--no-lammps           不运行LAMMPS模拟
--run-gaussian        运行高斯计算
--gaussian-method     高斯计算方法，默认为"B3LYP"
--gaussian-basis      高斯计算基组，默认为"6-31G(d)"
--gaussian-calc-type  高斯计算类型，默认为"opt freq"
```

### 输入Excel格式

输入Excel文件应遵循Molyte标准格式，包含以下必填字段：

1. **基本信息**：
   - `name`：系统名称
   - `temperature`：模拟温度 (K)
   - `box_size`：盒子尺寸
   - `concentration`：浓度

2. **组分信息**：
   - 阳离子：`cation1_name`、`cation1_smile`、`cation1_ratio`等
   - 阴离子：`anion1_name`、`anion1_smile`、`anion1_ratio`等
   - 溶剂：`sol1_name`、`sol1_smile`、`sol1_ratio`等

### 工作流程

1. 读取Excel文件，提取系统信息
2. 复制Excel文件到LAMMPS工作空间，添加时间戳
3. 对每个系统：
   - 设置并提交LAMMPS模拟作业
   - 如果需要，生成高斯输入文件并运行计算
4. 记录所有处理的系统信息

### 示例

```bash
# 基本用法
python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx

# 同时运行高斯计算
python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx --run-gaussian

# 只设置高斯计算，不运行LAMMPS
python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx --no-lammps --run-gaussian
```

## 后处理阶段

后处理阶段用于处理LAMMPS和高斯计算的结果，并整合到标准结果目录。

### 命令行参数

```
--postprocess          执行后处理阶段
--no-lammps-analysis   不运行LAMMPS结果分析
--no-gaussian-analysis 不运行分子计算结果分析
```

### 工作流程

1. 扫描工作目录结构，寻找所有LAMMPS和高斯计算结果
2. 对于LAMMPS结果：
   - 修复XYZ轨迹文件中的周期性边界问题
   - 复制关键结果文件到LAMMPS结果目录
3. 对于高斯计算结果：
   - 解析高斯输出文件，提取关键电子性质
   - 生成可视化图表
   - 复制结果文件到高斯结果目录

### 示例

```bash
# 处理所有结果
python -m molyte_cursor.src.core.integrated_workflow --postprocess

# 只处理LAMMPS结果
python -m molyte_cursor.src.core.integrated_workflow --postprocess --no-gaussian-analysis

# 只处理高斯计算结果
python -m molyte_cursor.src.core.integrated_workflow --postprocess --no-lammps-analysis
```

## 使用示例

### 完整工作流示例

以下是一个典型的完整工作流示例：

1. **准备输入Excel文件**
   - 创建包含系统信息的Excel文件
   - 确保所有必需字段都已填写

2. **执行预处理阶段**
   ```bash
   # 设置并运行LAMMPS模拟，同时生成高斯输入文件
   python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx --run-gaussian
   ```

3. **等待计算完成**
   - LAMMPS模拟可能需要几小时到几天
   - 高斯计算可能需要几分钟到几小时

4. **执行后处理阶段**
   ```bash
   # 分析所有结果
   python -m molyte_cursor.src.core.integrated_workflow --postprocess
   ```

5. **查看结果**
   - LAMMPS结果位于`Lammps_Results`目录
   - 高斯计算结果位于`Gaussian_Results`目录

### 只运行LAMMPS模拟

```bash
# 预处理阶段
python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx

# 后处理阶段
python -m molyte_cursor.src.core.integrated_workflow --postprocess --no-gaussian-analysis
```

### 只运行高斯计算

```bash
# 预处理阶段
python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx --no-lammps --run-gaussian

# 后处理阶段
python -m molyte_cursor.src.core.integrated_workflow --postprocess --no-lammps-analysis
```

## 常见问题

### 1. 预处理失败

问题: 预处理阶段失败，无法读取系统信息。

回答:
- 检查Excel文件格式是否正确
- 确保所有必需字段都已填写
- 查看日志输出以获取详细错误信息

### 2. LAMMPS模拟失败

问题: LAMMPS模拟设置或运行失败。

回答:
- 检查LAMMPS是否正确安装
- 确保系统参数在合理范围内
- 查看`Lammps_Workspace`中的日志文件获取详细错误信息

### 3. 高斯计算失败

问题: 高斯计算失败或未完成。

回答:
- 确保高斯软件正确安装
- 检查SMILE字符串是否有效
- 查看`Molecule_Workspace`中的日志文件获取详细错误信息

### 4. 结果文件未找到

问题: 后处理阶段找不到结果文件。

回答:
- 确保计算已完成
- 检查工作目录结构是否被手动修改
- 尝试手动指定工作空间目录：`--workspace /path/to/workspace`

### 5. 自定义工作空间

问题: 如何使用自定义工作空间？

回答:
- 使用`--workspace`参数指定工作空间根目录
- 例如：`python -m molyte_cursor.src.core.integrated_workflow --preprocess -i input.xlsx --workspace /path/to/workspace` 