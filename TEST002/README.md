# 材料属性计算与可视化平台

这是一个基于Flask和Pymatgen的材料属性计算与可视化平台，可以帮助研究人员快速分析晶体结构、计算材料属性，并生成VASP输入文件。

## 功能特点

- **结构可视化**：支持3D交互式结构查看（使用3Dmol.js和Three.js）
- **基本属性计算**：密度、晶胞体积、晶格参数、空间群等
- **高级属性分析**：
  - 最短原子距离和配位数计算
  - 表面能估算
  - 能带结构和态密度可视化（模拟数据）
- **VASP输入文件生成**：自动生成POSCAR、INCAR和KPOINTS文件

## 安装指南

### 环境要求

- Python 3.8+
- 依赖包：Flask, Pymatgen, NumPy等

### 安装步骤

1. 克隆仓库：
   ```bash
   git clone https://github.com/yourusername/material-properties-platform.git
   cd material-properties-platform
   ```

2. 创建虚拟环境（推荐）：
   ```bash
   python -m venv venv
   source venv/bin/activate  # Linux/Mac
   # 或
   venv\Scripts\activate  # Windows
   ```

3. 安装依赖：
   ```bash
   pip install -r requirements.txt
   ```

4. 启动应用：
   ```bash
   cd app
   python app.py
   ```

5. 在浏览器中访问：`http://localhost:5000`

## 使用说明

1. **材料属性计算**：
   - 上传结构文件（支持CIF、VASP POSCAR等格式）
   - 选择要计算的属性
   - 点击"计算属性"按钮

2. **结构可视化**：
   - 支持三种可视化模式：3Dmol、Three.js和表格视图
   - 可以旋转、缩放查看3D结构

3. **VASP输入文件生成**：
   - 选择"生成VASP输入文件"选项
   - 可以复制或下载生成的POSCAR、INCAR和KPOINTS文件

## 注意事项

- 表面能计算仅提供几何信息，不是真实的表面能（需要DFT计算）
- 能带结构和态密度使用模拟数据，不代表实际材料的电子结构
- 最短原子距离和配位数计算可能对某些复杂结构不准确

## 已知问题

- 在某些情况下，最短原子距离和配位数计算可能失败，这通常是由于结构复杂或Pymatgen的API变化导致的
- 表面能计算对某些结构可能不稳定

## 贡献指南

欢迎提交问题报告和改进建议！如果您想贡献代码，请遵循以下步骤：

1. Fork仓库
2. 创建功能分支：`git checkout -b feature/your-feature-name`
3. 提交更改：`git commit -m 'Add some feature'`
4. 推送到分支：`git push origin feature/your-feature-name`
5. 提交Pull Request

## 许可证

MIT License