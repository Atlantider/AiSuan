# Molyte轨迹处理工具使用指南

本文档介绍如何使用Molyte中的轨迹处理工具，特别是如何修复因周期性边界条件导致的分子断裂问题。

## 目录

- [问题背景](#问题背景)
- [解决方案](#解决方案)
- [命令行工具](#命令行工具)
- [编程接口](#编程接口)
- [常见问题](#常见问题)

## 问题背景

在分子动力学模拟中，通常使用周期性边界条件来避免表面效应并模拟宏观系统。然而，这种方法会在可视化或分析轨迹时产生一个问题：当分子跨越盒子边界时，在输出的坐标文件（如XYZ文件）中，这些分子会显得"断裂"，因为一部分原子出现在盒子的一侧，而另一部分出现在盒子的对侧。

这对溶剂化结构的分析尤其不利，因为溶剂分子应该紧密围绕溶质，但在可视化时看起来却是分散的。

## 解决方案

Molyte轨迹处理工具提供了解决这一问题的方法：

1. **识别分子**: 根据原子间距离确定哪些原子属于同一分子
2. **展开分子**: 将跨越周期性边界的分子重新连接在一起
3. **整体移动**: 确保每个分子都被正确地表示，不会因周期性边界而断裂

这个过程称为"unwrapping"（展开），是分子动力学后处理中的常见操作。

## 命令行工具

### 基本用法

```bash
# 处理单个XYZ文件
python -m molyte_cursor.src.analysis.fix_xyz path/to/trajectory.xyz

# 处理目录中的所有XYZ文件
python -m molyte_cursor.src.analysis.fix_xyz path/to/trajectory_directory/
```

### 命令行选项

```
-o, --output       指定输出文件路径（仅处理单个文件时有效）
-b, --bond-threshold  设置原子间连接的距离阈值（埃），默认为2.0
-w, --wrap-only    仅将原子包装回盒子内，不进行分子修复
-p, --pattern      处理目录时的文件匹配模式，默认为'*.xyz'
-v, --verbose      输出详细日志信息
```

### 使用示例

1. 处理单个轨迹文件并指定输出:
   ```bash
   python -m molyte_cursor.src.analysis.fix_xyz simulation.xyz -o fixed_simulation.xyz
   ```

2. 使用较大的键长阈值（对于某些带电体系可能需要):
   ```bash
   python -m molyte_cursor.src.analysis.fix_xyz simulation.xyz -b 2.5
   ```

3. 处理目录中所有XYZ文件:
   ```bash
   python -m molyte_cursor.src.analysis.fix_xyz trajectory_dir/ -v
   ```

4. 仅处理特定帧:
   ```bash
   python -m molyte_cursor.src.analysis.fix_xyz trajectory_dir/ -p "*frame_100.xyz"
   ```

## 编程接口

如果需要在Python脚本中使用该功能，可以直接导入`TrajectoryProcessor`类:

```python
from molyte_cursor.src.analysis.trajectory_processor import TrajectoryProcessor

# 创建轨迹处理器
processor = TrajectoryProcessor()

# 修复单个文件
processor.fix_periodic_molecules(
    "input.xyz",
    "output.xyz",
    bond_threshold=2.0,
    unwrap=True
)

# 处理整个目录
processor.process_trajectory_directory(
    "trajectory_dir/",
    file_pattern="*.xyz",
    unwrap=True
)
```

### 关键方法

- `fix_periodic_molecules`: 修复单个XYZ文件中的周期性问题
- `process_trajectory_directory`: 处理目录中的多个轨迹文件
- `unwrap_molecule`: 展开单个分子的坐标
- `wrap_coordinates`: 将坐标包装回模拟盒子内

## 常见问题

### 1. 如何确定合适的bond_threshold值?

问题: 默认的键长阈值(2.0埃)可能不适用于某些特殊体系，如何确定合适的值?

回答: 通常，2.0埃适用于大多数共价键。但对于某些体系:
- 对于含有过渡金属的体系，可能需要更大的值(2.2-2.5埃)
- 对于氢键或其他弱相互作用，可能需要更小的值(1.5-1.8埃)
- 通过观察修复前后的结构，可以调整这个值以获得最佳效果

### 2. 处理大型轨迹时内存问题

问题: 处理包含成千上万原子的轨迹文件时可能会消耗大量内存。

回答: 对于非常大的体系，可以:
- 先只处理几个关键帧进行测试
- 将轨迹分割成多个小文件分别处理
- 增加物理内存或使用配备更多RAM的机器

### 3. 无法识别某些分子

问题: 有时工具无法正确识别某些分子，特别是对于复杂的离子或大分子。

回答: 这通常是由于键长阈值设置不当。尝试:
- 调整`bond_threshold`参数
- 对于特殊体系，可能需要编写自定义的分子识别规则
- 检查输入XYZ文件中的原子坐标是否合理 