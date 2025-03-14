# Molyte Cursor

这是 Molyte 分子动力学模拟工具的重构版本，采用更规范的代码结构和组织方式。

## 项目结构

```
molyte_cursor/
├── config/                     # 配置文件目录
│   ├── paths/                  # 路径配置
│   │   └── default_paths.yaml  # 默认路径配置
│   └── simulation_params/      # 模拟参数配置
│       └── default_params.yaml # 默认模拟参数
├── src/                        # 源代码目录
│   ├── core/                   # 核心功能
│   │   └── main.py             # 主程序入口
│   ├── io/                     # 输入输出处理
│   │   ├── excel_reader.py     # Excel读取器
│   │   └── file_generator.py   # 文件生成器
│   ├── models/                 # 数据模型
│   │   └── molecule.py         # 分子模型类
│   ├── simulation/             # 模拟相关
│   │   └── simulator.py        # 模拟处理器
│   ├── analysis/               # 分析模块
│   │   ├── analyzer.py         # 分析主类
│   │   ├── rdf_analyzer.py     # RDF分析器
│   │   ├── msd_analyzer.py     # MSD分析器
│   │   ├── solvent_analyzer.py # 溶剂分析器
│   │   ├── calc_properties.py  # 物理化学性质计算
│   │   ├── visualization.py    # 可视化组件
│   │   └── README.md           # 分析模块使用指南
│   └── utils/                  # 工具类
│       ├── charge_modifier.py  # 电荷修改工具
│       ├── command_executor.py # 命令执行工具
│       ├── config_loader.py    # 配置加载工具
│       └── logger.py           # 日志工具
├── tests/                      # 测试目录
├── docs/                       # 文档目录
└── setup.py                    # 安装配置
```

## 功能特点

1. **模块化设计**: 将代码按功能模块分离，提高可维护性
2. **配置文件化**: 将固定参数和路径移至配置文件，便于修改
3. **规范化日志**: 统一的日志记录机制
4. **完整错误处理**: 强化错误处理和异常捕获
5. **面向对象设计**: 采用类和对象组织代码，提高可读性和复用性
6. **全面分析功能**: 提供RDF、MSD、溶剂分布、电化学性能等全面分析

## 安装方法

```bash
# 从源代码安装
pip install -e .
```

## 使用方法

```bash
# 基本用法
molyte --input /path/to/input.xlsx

# 提交计算作业
molyte --input /path/to/input.xlsx --submit

# 指定日志文件
molyte --input /path/to/input.xlsx --log /path/to/log.txt

# 指定配置目录
molyte --input /path/to/input.xlsx --config /path/to/config/dir
```

## 依赖工具

- LAMMPS
- Packmol
- Moltemplate
- ligpargen
- RESP

## 新增功能

### 电化学性能分析

最新版本新增了全面的电化学性能分析功能：

1. **离子电导率计算**：基于Nernst-Einstein方程计算电导率和摩尔电导率
2. **迁移数分析**：计算阳离子迁移数，评估电池电解液性能
3. **离子团簇分析**：识别接触离子对(CIP)、溶剂分离离子对(SSIP)和聚集体(AGG)
4. **粘度估算**：基于Stokes-Einstein关系从扩散系数估算粘度
5. **电化学数据可视化**：生成电导率、粘度、迁移数等专业图表
6. **多样品比较分析**：支持不同电解液配方在相同条件下的性能对比
7. **数据导出**：自动生成Excel报告，便于数据整理和分享

详细使用方法请参考 `src/analysis/README.md`

## 版本历史

### v1.2.0 (2023年3月)
- 添加电化学性能分析模块
- 新增离子团簇分析功能
- 改进批量分析和可视化组件

### v1.1.0 
- 优化代码结构，采用面向对象设计
- 添加配置文件管理系统

### v1.0.0
- 初始版本发布 