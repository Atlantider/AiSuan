# Molyte分析模块

本模块提供了分子动力学模拟结果的分析功能，包括RDF分析、MSD分析、高斯计算分析和溶剂分析等。通过面向对象的设计，使分析功能更加模块化和易于使用。

## 模块结构

- `RDFAnalyzer`: 径向分布函数分析器
- `MSDAnalyzer`: 均方位移分析器，包含扩散系数计算
- `GaussianAnalyzer`: 高斯量子化学计算分析器
- `SolventAnalyzer`: 溶剂分布和溶剂化层分析器
- `Visualizer`: 可视化器，支持各种分析结果的图形展示
- `Analyzer`: 分析运行器，整合所有分析器功能的主入口

## 使用方法

### 基本使用

```python
from molyte_cursor.src.analysis import Analyzer

# 初始化分析器
analyzer = Analyzer(
    base_dir='/path/to/simulation/results',  # 模拟结果所在的基础目录
    output_dir='/path/to/output'             # 分析结果输出目录
)

# 运行单个样品的RDF分析
rdf_results = analyzer.run_rdf_analysis('sample_name')

# 运行单个样品的MSD分析
msd_results = analyzer.run_msd_analysis('sample_name')

# 运行单个样品的溶剂分析
solvent_results = analyzer.run_solvent_analysis('sample_name')

# 运行综合分析（包括RDF、MSD和溶剂分析）
results = analyzer.run_comprehensive_analysis('sample_name')
```

### 批量分析

```python
# 批量分析多个样品
sample_names = ['sample1', 'sample2', 'sample3']
batch_results = analyzer.run_batch_analysis(
    sample_names=sample_names,
    analysis_types=['rdf', 'msd', 'solvent', 'comprehensive']
)

# 比较多个样品的MSD数据
comparative_results = analyzer.run_comparative_msd_analysis(sample_names)
```

### 高斯计算分析

```python
# 分析高斯计算结果
gaussian_results = analyzer.run_gaussian_analysis(
    directory='/path/to/gaussian/outputs',
    output_name='gaussian_analysis'
)
```

## 各分析器的单独使用

除了使用`Analyzer`整合类，你还可以单独使用各个分析器：

### RDF分析器

```python
from molyte_cursor.src.analysis import RDFAnalyzer

rdf_analyzer = RDFAnalyzer(working_dir='/path/to/base/dir')

# 解析in_list文件并提取RDF标签
element_lists, rdf_labels = rdf_analyzer.parse_in_list_with_molecule(
    in_list_path='sample.in.list',
    atom_counts={'cation': 20, 'anion': 30, 'solvent': 10}
)

# 替换RDF文件中的标签
rdf_analyzer.replace_rdf_labels_in_file(
    original_file='out_rdf.dat',
    rdf_labels=rdf_labels,
    output_file='rdf_parsed.dat'
)

# 绘制RDF图像
rdf_analyzer.plot_rdf(
    rdf_file='rdf_parsed.dat',
    rdf_labels=rdf_labels,
    output_filename='rdf_plot.png'
)
```

### MSD分析器

```python
from molyte_cursor.src.analysis import MSDAnalyzer

msd_analyzer = MSDAnalyzer(working_dir='/path/to/base/dir')

# 绘制MSD图像
msd_analyzer.plot_msd(
    msd_file='out_msd.dat',
    output_filename='msd_plot.png'
)

# 计算扩散系数
diffusion_coeffs = msd_analyzer.calculate_diffusion_coefficient(
    msd_file='out_msd.dat',
    time_unit='fs',
    skip_initial=10
)

# 处理单个样品的MSD数据
msd_analyzer.process_msd_data(
    sample_name='sample_name',
    current_directory='/path/to/sample/dir',
    output_directory='/path/to/output/dir'
)
```

### 溶剂分析器

```python
from molyte_cursor.src.analysis import SolventAnalyzer

solvent_analyzer = SolventAnalyzer(working_dir='/path/to/base/dir')

# 分析溶剂分布
solvent_results = solvent_analyzer.analyze_solvent_distribution(
    dump_file='trajectory.dump',
    atom_counts={'cation': 20, 'anion': 30, 'solvent': 10},
    cation_atom_types=[1, 2],
    anion_atom_types=[3, 4],
    solvent_atom_types=[5, 6],
    cutoff_distance=5.0
)
```

### 可视化器

```python
from molyte_cursor.src.analysis import Visualizer

visualizer = Visualizer(
    output_dir='/path/to/visualizations',
    dpi=300
)

# 绘制RDF图像
visualizer.plot_rdf(
    rdf_file='rdf_parsed.dat',
    rdf_labels=rdf_labels
)

# 绘制MSD图像
visualizer.plot_msd(
    msd_file='out_msd.dat'
)

# 绘制扩散系数对比图
visualizer.plot_diffusion_coefficients(
    diffusion_data={
        'sample1': {'D_x': 1.2, 'D_y': 1.3, 'D_z': 1.4, 'D_total': 1.3},
        'sample2': {'D_x': 2.2, 'D_y': 2.3, 'D_z': 2.4, 'D_total': 2.3}
    },
    output_filename='diffusion_comparison.png'
)

# 绘制溶剂分布图
visualizer.plot_solvent_distribution(
    solvent_data=solvent_results,
    output_filename='solvent_distribution.png',
    plot_type='bar'  # 可选：'bar', 'pie', 'heatmap'
)
```

## 输出结果

所有分析结果都将保存在指定的输出目录中，包括：

1. 分析数据文件（JSON格式）
2. 图像文件（PNG格式）
3. 处理后的数据文件（DAT、TXT格式）
4. 表格数据（XLSX格式）

每个样品的分析结果将保存在以样品名称命名的子目录中，便于查找和管理。

## 电化学性能分析

电化学性能分析功能可以基于MSD数据计算离子电导率、迁移数和粘度等电化学参数，为研究电解液性能提供重要参考。

### 单样品电化学分析

以下示例展示了如何对单个样品进行电化学性能分析：

```python
from molyte_cursor.src.analysis import Analyzer

# 初始化分析器
analyzer = Analyzer(
    base_dir="/path/to/data_directory",
    output_dir="/path/to/output_directory"
)

# 运行电化学性能分析
# 注意：这需要先运行MSD分析以获取必要的MSD数据
results = analyzer.run_electrochemical_analysis(
    sample_name="LiTFSI-EC",
    temperature=298.15,  # 温度 (K)
    concentration=1.0,   # 电解质浓度 (mol/L)
    cation_valence=1,    # 阳离子价数 (Li+ = 1)
    anion_valence=1,     # 阴离子价数 (TFSI- = 1)
    particle_radius=2.0  # 粒子半径 (Å)，用于估算粘度
)

# 结果包含以下参数：
# - conductivity_mS_cm: 电导率 (mS/cm)
# - molar_conductivity: 摩尔电导率 (S·m²/mol)
# - transference_number: 阳离子迁移数
# - viscosity_mPa_s: 粘度 (mPa·s)
# - diffusion_coeffs: 各组分扩散系数 (10⁻⁹ m²/s)
```

### 多样品对比分析

以下示例展示了如何对多个样品的电化学性能进行对比分析：

```python
from molyte_cursor.src.analysis import Analyzer

# 初始化分析器
analyzer = Analyzer(
    base_dir="/path/to/data_directory",
    output_dir="/path/to/output_directory"
)

# 运行多样品对比分析
results = analyzer.run_comparative_electrochemical_analysis(
    sample_names=["LiTFSI-EC", "LiTFSI-DMC", "LiTFSI-EC-DMC"],
    output_name="lithium_electrolytes_comparison",
    temperature=298.15,  # 温度 (K)
    concentrations=[1.0, 1.0, 1.0],  # 各样品的电解质浓度 (mol/L)
    cation_valence=1,    # 阳离子价数
    anion_valence=1      # 阴离子价数
)

# 结果将保存为Excel文件，并生成各项电化学参数的对比图表
```

### 批量分析

以下示例展示了如何在批量分析中包含电化学性能分析：

```python
from molyte_cursor.src.analysis import Analyzer

# 初始化分析器
analyzer = Analyzer(
    base_dir="/path/to/data_directory",
    output_dir="/path/to/output_directory"
)

# 批量运行多种分析，包括电化学分析
results = analyzer.run_batch_analysis(
    sample_names=["LiTFSI-EC", "LiTFSI-DMC", "LiTFSI-EC-DMC"],
    analysis_types=["msd", "electrochemical"]  # 先运行MSD分析，再进行电化学分析
)
```

### 主要电化学参数说明

1. **离子电导率**：表示电解液导电能力的关键参数，单位为mS/cm。通过Nernst-Einstein方程从离子扩散系数计算得到。

2. **摩尔电导率**：考虑了浓度因素的电导率标准化值，单位为S·m²/mol。摩尔电导率与浓度的关系可用于研究离子聚集现象。

3. **阳离子迁移数**：表示电流中由阳离子迁移贡献的比例，范围为0-1。理想的锂离子电池电解液应有接近0.5或更高的锂离子迁移数。

4. **粘度**：表示流体阻力的物理量，单位为mPa·s。粘度与扩散系数成反比，通过Stokes-Einstein关系从扩散系数估算。 

## 高级电化学计算函数

除了通过主分析器使用电化学分析功能外，您还可以单独使用各种电化学性能计算函数：

### 电导率计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_ionic_conductivity

# 从MSD数据计算离子电导率
conductivity = calculate_ionic_conductivity(
    msd_file="path/to/out_msd.dat",          # MSD数据文件路径
    temperature=298.15,                       # 温度 (K)
    time_unit='fs',                           # 时间单位
    concentration=1.0,                        # 浓度 (mol/L)
    cation_valence=1,                         # 阳离子价数
    anion_valence=1,                          # 阴离子价数
    skip_initial=10                           # 忽略初始数据点数量
)

# 返回值单位为S/m，可转换为mS/cm（乘以10）
```

### 介电常数计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_dielectric_constant

# 计算介电常数
dielectric = calculate_dielectric_constant(
    dipole_file="path/to/dipole.dat",         # 偶极矩数据文件路径
    temperature=298.15,                       # 温度 (K)
    volume=1000.0,                            # 系统体积 (Å³)
    dipole_unit='e*A'                         # 偶极矩单位
)
```

### 迁移数计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_transference_number

# 计算阳离子迁移数
t_plus = calculate_transference_number(
    cation_diffusion=1.2e-9,                  # 阳离子扩散系数 (m²/s)
    anion_diffusion=0.8e-9,                   # 阴离子扩散系数 (m²/s)
    cation_valence=1,                         # 阳离子价数
    anion_valence=1                           # 阴离子价数
)
```

### 哈文比（Haven Ratio）计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_haven_ratio

# 计算哈文比（离子关联度）
HR = calculate_haven_ratio(
    conductivity_md=0.5,                      # MD模拟计算的电导率
    conductivity_ne=0.8                       # 根据Nernst-Einstein方程计算的电导率
)
```

### 扩散激活能计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_diffusion_activation_energy

# 计算扩散激活能
Ea, A = calculate_diffusion_activation_energy(
    diffusion_coeffs=[1.0e-9, 1.5e-9, 2.5e-9],    # 不同温度下的扩散系数列表
    temperatures=[273.15, 298.15, 323.15]         # 对应的温度列表 (K)
)
# Ea为激活能（kJ/mol），A为指前因子
```

### 粘度计算

```python
from molyte_cursor.src.analysis.calc_properties import calculate_viscosity

# 从扩散系数估算粘度
viscosity = calculate_viscosity(
    diffusion_coeff=1.2e-9,                   # 扩散系数 (m²/s)
    particle_radius=2.0,                      # 粒子半径 (Å)
    temperature=298.15                        # 温度 (K)
)
# 返回值单位为Pa·s，可转换为mPa·s（乘以1000）
```

## 电化学数据可视化

我们提供了一系列函数用于可视化电化学分析结果，可以通过Visualizer类直接调用：

```python
from molyte_cursor.src.analysis import Visualizer

visualizer = Visualizer(output_dir='/path/to/output')

# 绘制单个样品的电导率数据
conductivity_plot = visualizer.plot_conductivity_data(
    electrochemical_data={
        'sample_name': 'LiTFSI-EC',
        'conductivity_mS_cm': 5.43,
        'molar_conductivity': 0.00543,
        'concentration': 1.0
    },
    output_filename='conductivity_plot.png'
)

# 绘制单个样品的粘度数据
viscosity_plot = visualizer.plot_viscosity_data(
    electrochemical_data={
        'sample_name': 'LiTFSI-EC',
        'viscosity_mPa_s': 3.21,
        'diffusion_coeffs': {'D_total': 1.25}
    },
    output_filename='viscosity_plot.png'
)

# 绘制多个样品的电导率比较图
comparative_plot = visualizer.plot_comparative_conductivity(
    electrochemical_data={
        'LiTFSI-EC': {
            'conductivity_mS_cm': 5.43,
            'molar_conductivity': 0.00543
        },
        'LiTFSI-DMC': {
            'conductivity_mS_cm': 7.65,
            'molar_conductivity': 0.00765
        }
    },
    output_filename='comparative_conductivity.png'
)

# 绘制多个样品的粘度比较图
viscosity_plot = visualizer.plot_comparative_viscosity(
    electrochemical_data={
        'LiTFSI-EC': {
            'viscosity_mPa_s': 3.21,
            'diffusion_coefficients': {'D_total': 1.25}
        },
        'LiTFSI-DMC': {
            'viscosity_mPa_s': 2.14,
            'diffusion_coefficients': {'D_total': 1.87}
        }
    },
    output_filename='comparative_viscosity.png'
)

# 绘制多个样品的迁移数比较图
transference_plot = visualizer.plot_comparative_transference(
    electrochemical_data={
        'LiTFSI-EC': {'transference_number': 0.35},
        'LiTFSI-DMC': {'transference_number': 0.42},
        'LiTFSI-EC-DMC': {'transference_number': 0.38}
    },
    output_filename='comparative_transference.png'
)
```

## 电化学分析输出数据格式

电化学分析完成后，会生成多种格式的输出文件：

### JSON结果文件

包含完整分析结果的JSON文件，包含以下主要键值：

```json
{
  "sample_name": "样品名称",
  "temperature": 298.15,
  "concentration": 1.0,
  "diffusion_coeffs": {
    "D_total": 1.25,
    "D_cation": 1.0,
    "D_anion": 0.8,
    "D_solvent": 1.5,
    "D_x": 1.2,
    "D_y": 1.3,
    "D_z": 1.25
  },
  "conductivity": 0.543,
  "conductivity_mS_cm": 5.43,
  "molar_conductivity": 0.00543,
  "transference_number": 0.35,
  "viscosity": 0.00321,
  "viscosity_mPa_s": 3.21,
  "plots": {
    "conductivity_plot": "path/to/conductivity_plot.png",
    "viscosity_plot": "path/to/viscosity_plot.png"
  }
}
```

### Excel结果文件

比较分析会生成Excel文件，包含以下工作表：

1. **比较结果**: 包含多个样品的对比数据表格
2. **图表列表**: 生成的图表文件路径

### 图像文件

生成的图像包括：

1. **电导率图**: 显示电导率和摩尔电导率
2. **粘度图**: 显示粘度和总扩散系数
3. **迁移数图**: 显示阳离子迁移数

所有图像均为高分辨率PNG格式，默认DPI为300。 