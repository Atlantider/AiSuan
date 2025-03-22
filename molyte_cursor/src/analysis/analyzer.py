"""
分析运行器模块

提供了运行各种分析任务的功能，整合了所有分析器的功能
"""

import os
import logging
import json
from typing import Dict, List, Tuple, Optional, Union, Any
from collections import defaultdict
import pandas as pd

from .rdf_analyzer import RDFAnalyzer
from .msd_analyzer import MSDAnalyzer
from .gaussian_analyzer import GaussianAnalyzer
from .solvent_analyzer import SolventAnalyzer
from .visualization import Visualizer

class Analyzer:
    """
    分析运行器
    
    整合了所有分析器的功能，提供了运行各种分析任务的接口
    """
    def __init__(self, base_dir: str, output_dir: Optional[str] = None, logger=None):
        """
        初始化分析运行器

        Args:
            base_dir: 基础目录，存放模拟结果
            output_dir: 输出目录，存放分析结果，如果为None则使用base_dir/analysis
            logger: 日志记录器实例
        """
        self.base_dir = os.path.abspath(base_dir)
        self.output_dir = os.path.abspath(output_dir) if output_dir else os.path.join(self.base_dir, 'analysis')
        self.logger = logger or logging.getLogger(__name__)
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 初始化各个分析器
        self.rdf_analyzer = RDFAnalyzer(self.base_dir, self.logger)
        self.msd_analyzer = MSDAnalyzer(self.base_dir, self.logger)
        self.gaussian_analyzer = GaussianAnalyzer(self.base_dir, 
                                                 os.path.join(self.output_dir, 'energy_cache.json'),
                                                 self.logger)
        self.solvent_analyzer = SolventAnalyzer()
        self.visualizer = Visualizer(os.path.join(self.output_dir, 'visualizations'), 
                                    dpi=300, logger=self.logger)
    
    def run_rdf_analysis(self, sample_name: str, in_list_path: Optional[str] = None, 
                        atom_counts: Optional[Dict[str, int]] = None) -> Dict[str, Any]:
        """
        运行RDF分析

        Args:
            sample_name: 样品名称
            in_list_path: 输入列表文件路径，如果为None则使用默认路径
            atom_counts: 分子原子数量映射，如果为None则尝试加载

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行RDF分析: {sample_name}")
        
        # 确定目录和文件路径
        sample_dir = os.path.join(self.base_dir, sample_name)
        output_sample_dir = os.path.join(self.output_dir, sample_name)
        os.makedirs(output_sample_dir, exist_ok=True)
        
        # 设置默认路径
        if in_list_path is None:
            in_list_path = os.path.join(sample_dir, f"{sample_name}.in.list")
        
        # 如果未提供atom_counts，尝试从文件加载
        if atom_counts is None:
            atom_counts_file = os.path.join(sample_dir, 'atom_counts.json')
            if os.path.exists(atom_counts_file):
                with open(atom_counts_file, 'r') as f:
                    atom_counts = json.load(f)
            else:
                self.logger.error(f"未提供atom_counts且找不到atom_counts文件: {atom_counts_file}")
                return {'error': 'Missing atom_counts'}
        
        # 原始RDF文件路径
        original_rdf_file = os.path.join(sample_dir, 'out_rdf.dat')
        if not os.path.exists(original_rdf_file):
            self.logger.error(f"找不到RDF文件: {original_rdf_file}")
            return {'error': 'Missing RDF file'}
        
        # 解析in_list文件并提取元素列表和RDF标签
        element_lists, rdf_labels = self.rdf_analyzer.parse_in_list_with_molecule(in_list_path, atom_counts)
        
        # 替换RDF标签
        rdf_output_path = os.path.join(output_sample_dir, f'rdf_{sample_name}_parser.dat')
        self.rdf_analyzer.replace_rdf_labels_in_file(original_rdf_file, rdf_labels, rdf_output_path)
        
        # 绘制RDF图像
        rdf_image_path = os.path.join(output_sample_dir, f'rdf_{sample_name}.png')
        self.visualizer.plot_rdf(rdf_output_path, rdf_labels, rdf_image_path)
        
        results = {
            'sample_name': sample_name,
            'element_lists': element_lists,
            'rdf_labels': rdf_labels,
            'rdf_output_path': rdf_output_path,
            'rdf_image_path': rdf_image_path
        }
        
        # 保存结果到JSON文件
        result_json_path = os.path.join(output_sample_dir, f'rdf_analysis_{sample_name}.json')
        with open(result_json_path, 'w') as f:
            # 过滤掉不可序列化的对象
            serializable_results = {k: v for k, v in results.items() 
                                   if isinstance(v, (str, int, float, list, dict, bool, type(None)))}
            json.dump(serializable_results, f, indent=2)
        
        self.logger.info(f"RDF分析完成: {sample_name}")
        return results
    
    def run_msd_analysis(self, sample_name: str) -> Dict[str, Any]:
        """
        运行MSD分析

        Args:
            sample_name: 样品名称

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行MSD分析: {sample_name}")
        
        # 确定目录和文件路径
        sample_dir = os.path.join(self.base_dir, sample_name)
        output_sample_dir = os.path.join(self.output_dir, sample_name)
        os.makedirs(output_sample_dir, exist_ok=True)
        
        # 处理MSD数据
        diffusion_coeffs = self.msd_analyzer.process_msd_data(sample_name, sample_dir, output_sample_dir)
        
        if not diffusion_coeffs:
            self.logger.warning(f"未找到MSD数据或处理失败: {sample_name}")
            return {'error': 'MSD data processing failed'}
        
        results = {
            'sample_name': sample_name,
            'diffusion_coefficients': diffusion_coeffs,
            'msd_image_path': os.path.join(output_sample_dir, f'msd_{sample_name}.png'),
            'diffusion_file': os.path.join(output_sample_dir, f'diffusion_{sample_name}.txt')
        }
        
        # 保存结果到JSON文件
        result_json_path = os.path.join(output_sample_dir, f'msd_analysis_{sample_name}.json')
        with open(result_json_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        self.logger.info(f"MSD分析完成: {sample_name}")
        return results
    
    def run_comparative_msd_analysis(self, sample_names: List[str], output_name: str = 'combined_msd') -> Dict[str, Any]:
        """
        运行比较性MSD分析，对多个样品进行比较

        Args:
            sample_names: 样品名称列表
            output_name: 输出文件名前缀

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行比较性MSD分析: {sample_names}")
        
        # 收集所有样品的MSD文件和扩散系数
        msd_files = {}
        diffusion_data = {}
        
        for sample_name in sample_names:
            sample_dir = os.path.join(self.base_dir, sample_name)
            msd_file = os.path.join(sample_dir, 'out_msd.dat')
            
            if os.path.exists(msd_file):
                msd_files[sample_name] = msd_file
                
                # 尝试加载或计算扩散系数
                output_sample_dir = os.path.join(self.output_dir, sample_name)
                diffusion_file = os.path.join(output_sample_dir, f'diffusion_{sample_name}.txt')
                
                if os.path.exists(diffusion_file):
                    # 从文件加载扩散系数
                    try:
                        with open(diffusion_file, 'r') as f:
                            lines = f.readlines()
                            diffusion_coeffs = {}
                            for line in lines[2:]:  # 跳过前两行
                                if ':' in line:
                                    key, value = line.strip().split(':')
                                    key = key.strip()
                                    value = float(value.strip())
                                    diffusion_coeffs[f"D_{key.split('方向')[0].lower()}"] = value
                            diffusion_data[sample_name] = diffusion_coeffs
                    except Exception as e:
                        self.logger.warning(f"加载扩散系数文件失败: {diffusion_file}, {e}")
                        # 如果加载失败，尝试重新计算
                        diffusion_coeffs = self.msd_analyzer.calculate_diffusion_coefficient(msd_file)
                        diffusion_data[sample_name] = diffusion_coeffs
                else:
                    # 计算扩散系数
                    os.makedirs(output_sample_dir, exist_ok=True)
                    diffusion_coeffs = self.msd_analyzer.calculate_diffusion_coefficient(msd_file)
                    diffusion_data[sample_name] = diffusion_coeffs
                    
                    # 保存扩散系数到文件
                    with open(diffusion_file, 'w') as f:
                        f.write(f"样品: {sample_name}\n")
                        f.write("扩散系数 (10^-9 m^2/s):\n")
                        f.write(f"X方向: {diffusion_coeffs['D_x']:.6f}\n")
                        f.write(f"Y方向: {diffusion_coeffs['D_y']:.6f}\n")
                        f.write(f"Z方向: {diffusion_coeffs['D_z']:.6f}\n")
                        f.write(f"总体: {diffusion_coeffs['D_total']:.6f}\n")
            else:
                self.logger.warning(f"找不到MSD文件: {msd_file}")
        
        if not msd_files:
            self.logger.error("没有找到任何有效的MSD文件")
            return {'error': 'No valid MSD files found'}
        
        # 生成组合MSD图
        combined_msd_img = self.visualizer.plot_combined_msd(
            msd_files, 
            f"{output_name}.png"
        )
        
        # 生成扩散系数对比图
        diffusion_img = self.visualizer.plot_diffusion_coefficients(
            diffusion_data, 
            f"{output_name}_diffusion.png"
        )
        
        results = {
            'sample_names': sample_names,
            'msd_files': msd_files,
            'diffusion_data': diffusion_data,
            'combined_msd_image': combined_msd_img,
            'diffusion_image': diffusion_img
        }
        
        # 保存结果到JSON文件
        result_json_path = os.path.join(self.output_dir, f'{output_name}_analysis.json')
        with open(result_json_path, 'w') as f:
            # 过滤掉不可序列化的对象
            serializable_results = {k: v for k, v in results.items() 
                                   if isinstance(v, (str, int, float, list, dict, bool, type(None)))}
            json.dump(serializable_results, f, indent=2)
        
        self.logger.info(f"比较性MSD分析完成")
        return results
    
    def run_electrochemical_analysis(self, sample_name: str, temperature: float, concentration: float,
                                   cation_valence: int = 1, anion_valence: int = 1,
                                   particle_radius: float = 2.0) -> Dict[str, Any]:
        """
        运行电化学性能分析，包括离子电导率、迁移数和粘度分析。

        Args:
            sample_name: 样品名称
            temperature: 温度 (K)
            concentration: 电解质浓度 (mol/L)
            cation_valence: 阳离子价数
            anion_valence: 阴离子价数
            particle_radius: 粒子半径 (Å) 用于估算粘度

        Returns:
            结果字典，包含电导率、迁移数、粘度等信息
        """
        self.logger.info(f"运行电化学性能分析: 样品名称={sample_name}, 温度={temperature}K")
        
        # 创建输出目录
        output_dir = os.path.join(self.output_dir, sample_name, "electrochemical_analysis")
        os.makedirs(output_dir, exist_ok=True)
        
        # 处理电化学数据
        msd_output_path = os.path.join(self.output_dir, sample_name, "msd_analysis", "out_msd.dat")
        if not os.path.exists(msd_output_path):
            self.logger.warning(f"未找到MSD数据文件: {msd_output_path}")
            self.logger.warning("请先运行MSD分析以生成所需的数据文件")
            return {"error": "未找到MSD数据文件"}
        
        # 处理MSD数据并计算扩散系数、电导率等
        electrochemical_results = self.msd_analyzer.process_electrochemical_data(
            sample_name=sample_name,
            msd_file=msd_output_path,
            output_dir=output_dir,
            temperature=temperature,
            concentration=concentration,
            cation_valence=cation_valence,
            anion_valence=anion_valence,
            particle_radius=particle_radius
        )
        
        # 生成电导率图表
        conductivity_plot_path = self.visualizer.plot_conductivity_data(
            electrochemical_data=electrochemical_results,
            output_filename=f"{sample_name}_conductivity.png"
        )
        
        # 生成粘度图表
        viscosity_plot_path = self.visualizer.plot_viscosity_data(
            electrochemical_data=electrochemical_results,
            output_filename=f"{sample_name}_viscosity.png"
        )
        
        # 将图表路径添加到结果中
        electrochemical_results["plots"] = {
            "conductivity_plot": conductivity_plot_path,
            "viscosity_plot": viscosity_plot_path
        }
        
        # 将结果保存为JSON
        results_json_path = os.path.join(output_dir, f"{sample_name}_electrochemical_results.json")
        with open(results_json_path, 'w') as f:
            json.dump(electrochemical_results, f, indent=4)
        
        self.logger.info(f"电化学性能分析完成: {sample_name}")
        self.logger.info(f"结果已保存到: {results_json_path}")
        
        # 返回结果和路径
        return {
            **electrochemical_results,
            "output_dir": output_dir,
            "results_json": results_json_path
        }
    
    def run_comparative_electrochemical_analysis(self, sample_names: List[str], output_name: str,
                                                   temperature: float, concentrations: Union[float, List[float]] = None,
                                                   cation_valence: int = 1, anion_valence: int = 1) -> Dict[str, Any]:
        """
        运行比较性电化学性能分析，对多个样品进行比较。

        Args:
            sample_names: 要比较的样品名称列表
            output_name: 输出文件名称前缀
            temperature: 温度 (K)
            concentrations: 电解质浓度 (mol/L)，可以是单一值或列表（与样品一一对应）
            cation_valence: 阳离子价数
            anion_valence: 阴离子价数

        Returns:
            结果字典，包含比较分析结果和图表路径
        """
        self.logger.info(f"运行比较性电化学性能分析: 样品={sample_names}, 温度={temperature}K")
        
        # 创建输出目录
        output_dir = os.path.join(self.output_dir, "comparative_analysis", output_name)
        os.makedirs(output_dir, exist_ok=True)
        
        # 处理浓度参数
        if concentrations is None:
            concentrations = [1.0] * len(sample_names)
        elif isinstance(concentrations, (int, float)):
            concentrations = [float(concentrations)] * len(sample_names)
        
        if len(concentrations) != len(sample_names):
            self.logger.warning(f"浓度列表长度 ({len(concentrations)}) 与样品数量 ({len(sample_names)}) 不匹配，使用默认值1.0")
            concentrations = [1.0] * len(sample_names)
        
        # 收集每个样品的电化学数据
        results = {}
        all_electrochemical_data = {}
        
        for i, sample_name in enumerate(sample_names):
            concentration = concentrations[i]
            
            self.logger.info(f"处理样品 {sample_name} (浓度: {concentration} mol/L)")
            
            # 运行单个样品的电化学分析
            sample_results = self.run_electrochemical_analysis(
                sample_name=sample_name,
                temperature=temperature,
                concentration=concentration,
                cation_valence=cation_valence,
                anion_valence=anion_valence
            )
            
            # 如果分析成功，收集结果
            if "error" not in sample_results:
                results[sample_name] = sample_results
                all_electrochemical_data[sample_name] = {
                    "sample_name": sample_name,
                    "concentration": concentration,
                    "conductivity_mS_cm": sample_results.get("conductivity_mS_cm", 0),
                    "molar_conductivity": sample_results.get("molar_conductivity", 0),
                    "transference_number": sample_results.get("transference_number", 0),
                    "viscosity_mPa_s": sample_results.get("viscosity_mPa_s", 0),
                    "diffusion_coefficients": sample_results.get("diffusion_coeffs", {})
                }
            else:
                self.logger.warning(f"样品 {sample_name} 分析失败: {sample_results.get('error')}")
        
        # 如果没有成功分析任何样品，返回错误
        if not results:
            self.logger.error("没有成功分析任何样品")
            return {"error": "没有成功分析任何样品"}
        
        # 生成比较图表
        plots = {}
        
        # 电导率比较图表
        conductivity_plot_path = self.visualizer.plot_comparative_conductivity(
            electrochemical_data=all_electrochemical_data,
            output_filename=f"{output_name}_comparative_conductivity.png"
        )
        plots["comparative_conductivity"] = conductivity_plot_path
        
        # 粘度比较图表
        viscosity_plot_path = self.visualizer.plot_comparative_viscosity(
            electrochemical_data=all_electrochemical_data,
            output_filename=f"{output_name}_comparative_viscosity.png"
        )
        plots["comparative_viscosity"] = viscosity_plot_path
        
        # 离子迁移数比较图表
        transference_plot_path = self.visualizer.plot_comparative_transference(
            electrochemical_data=all_electrochemical_data,
            output_filename=f"{output_name}_comparative_transference.png"
        )
        plots["comparative_transference"] = transference_plot_path
        
        # 创建Excel文件汇总结果
        excel_path = os.path.join(output_dir, f"{output_name}_comparative_results.xlsx")
        
        try:
            # 创建DataFrame
            data = []
            for sample, result in all_electrochemical_data.items():
                row = {
                    "样品名称": sample,
                    "浓度 (mol/L)": result.get("concentration", 0),
                    "电导率 (mS/cm)": result.get("conductivity_mS_cm", 0),
                    "摩尔电导率 (S·m²/mol)": result.get("molar_conductivity", 0),
                    "阳离子迁移数": result.get("transference_number", 0),
                    "粘度 (mPa·s)": result.get("viscosity_mPa_s", 0),
                    "总扩散系数 (10⁻⁹ m²/s)": result.get("diffusion_coefficients", {}).get("D_total", 0),
                    "阳离子扩散系数 (10⁻⁹ m²/s)": result.get("diffusion_coefficients", {}).get("D_cation", 0),
                    "阴离子扩散系数 (10⁻⁹ m²/s)": result.get("diffusion_coefficients", {}).get("D_anion", 0),
                    "溶剂扩散系数 (10⁻⁹ m²/s)": result.get("diffusion_coefficients", {}).get("D_solvent", 0)
                }
                data.append(row)
            
            df = pd.DataFrame(data)
            
            # 将结果写入Excel文件
            with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                df.to_excel(writer, sheet_name='比较结果', index=False)
                
                # 添加图表描述页
                plot_desc = pd.DataFrame({
                    "图表名称": list(plots.keys()),
                    "图表路径": list(plots.values())
                })
                plot_desc.to_excel(writer, sheet_name='图表列表', index=False)
            
            self.logger.info(f"比较结果已保存到Excel文件: {excel_path}")
            
        except Exception as e:
            self.logger.error(f"创建Excel结果文件时出错: {e}")
        
        # 返回结果和路径
        comparative_results = {
            "samples": sample_names,
            "sample_results": results,
            "plots": plots,
            "excel_file": excel_path,
            "output_dir": output_dir
        }
        
        return comparative_results
    
    def run_gaussian_analysis(self, directory: str, output_name: str = 'gaussian_analysis') -> Dict[str, Any]:
        """
        运行高斯计算分析

        Args:
            directory: 包含高斯输出文件的目录
            output_name: 输出文件名前缀

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行高斯计算分析: {directory}")
        
        if not os.path.isdir(directory):
            directory = os.path.join(self.base_dir, directory)
            if not os.path.isdir(directory):
                self.logger.error(f"找不到高斯输出目录: {directory}")
                return {'error': 'Invalid directory'}
        
        output_dir = os.path.join(self.output_dir, 'gaussian')
        os.makedirs(output_dir, exist_ok=True)
        
        # 分析所有高斯输出
        analysis_results = self.gaussian_analyzer.analyze_all_gaussian_outputs(directory)
        
        if not analysis_results:
            self.logger.warning(f"在目录中未找到有效的高斯输出文件: {directory}")
            return {'error': 'No valid Gaussian output files found'}
        
        # 提取能量数据和热力学数据
        energy_data = {}
        thermo_data = {}
        
        for filename, data in analysis_results.items():
            if data['finished'] and data['energy'] is not None:
                molecule_name = data['basename']
                energy_data[molecule_name] = data['energy']
                
                # 提取热力学数据
                thermo_results = self.gaussian_analyzer.calculate_thermochemistry(data['path'])
                if thermo_results:
                    thermo_data[molecule_name] = thermo_results
        
        # 可视化能量数据
        energy_img = None
        if energy_data:
            energy_img = self.visualizer.plot_gaussian_energies(
                energy_data, 
                f"{output_name}_energies.png"
            )
        
        # 可视化热力学数据
        thermo_img = None
        if thermo_data:
            thermo_img = self.visualizer.plot_thermochemistry(
                thermo_data, 
                f"{output_name}_thermochemistry.png"
            )
        
        results = {
            'directory': directory,
            'analysis_results': analysis_results,
            'energy_data': energy_data,
            'thermo_data': thermo_data,
            'energy_image': energy_img,
            'thermo_image': thermo_img
        }
        
        # 保存结果到JSON文件
        result_json_path = os.path.join(output_dir, f'{output_name}.json')
        with open(result_json_path, 'w') as f:
            # 过滤掉不可序列化的对象
            serializable_results = {k: v for k, v in results.items() 
                                   if isinstance(v, (str, int, float, list, dict, bool, type(None)))}
            json.dump(serializable_results, f, indent=2)
        
        self.logger.info(f"高斯计算分析完成")
        return results
    
    def run_solvent_analysis(self, sample_name: str, dump_file: Optional[str] = None,
                           atom_counts: Optional[Dict[str, int]] = None,
                           cation_atom_types: Optional[List[int]] = None,
                           anion_atom_types: Optional[List[int]] = None,
                           solvent_atom_types: Optional[List[int]] = None,
                           cutoff_distance: float = 5.0,
                           cip_distance: float = 4.0,
                           ssip_distance: float = 7.0) -> Dict[str, Any]:
        """
        运行溶剂分析

        Args:
            sample_name: 样品名称
            dump_file: LAMMPS轨迹文件路径，如果为None则使用默认路径
            atom_counts: 分子原子数量映射，如果为None则尝试加载
            cation_atom_types: 阳离子中代表性原子的类型列表
            anion_atom_types: 阴离子中代表性原子的类型列表
            solvent_atom_types: 溶剂中代表性原子的类型列表
            cutoff_distance: 溶剂化层的截止距离（埃）
            cip_distance: CIP的最大距离阈值（埃）
            ssip_distance: SSIP的最大距离阈值（埃）

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行溶剂分析: {sample_name}")
        
        # 确定目录和文件路径
        sample_dir = os.path.join(self.base_dir, sample_name)
        output_sample_dir = os.path.join(self.output_dir, sample_name)
        os.makedirs(output_sample_dir, exist_ok=True)
        
        # 设置默认路径
        if dump_file is None:
            dump_file = os.path.join(sample_dir, f"{sample_name}.dump")
            if not os.path.exists(dump_file):
                dump_file = os.path.join(sample_dir, 'dump.lammpstrj')
        
        if not os.path.exists(dump_file):
            self.logger.error(f"找不到LAMMPS轨迹文件: {dump_file}")
            return {'error': 'Missing LAMMPS trajectory file'}
        
        # 如果未提供atom_counts，尝试从文件加载
        if atom_counts is None:
            atom_counts_file = os.path.join(sample_dir, 'atom_counts.json')
            if os.path.exists(atom_counts_file):
                with open(atom_counts_file, 'r') as f:
                    atom_counts = json.load(f)
            else:
                self.logger.error(f"未提供atom_counts且找不到atom_counts文件: {atom_counts_file}")
                return {'error': 'Missing atom_counts'}
        
        # 尝试自动检测原子类型
        if cation_atom_types is None or anion_atom_types is None or solvent_atom_types is None:
            # 这里可以实现自动检测的逻辑
            # 作为简单示例，我们使用一些默认值
            cation_atom_types = cation_atom_types or [1, 2]  # 默认阳离子原子类型
            anion_atom_types = anion_atom_types or [3, 4]    # 默认阴离子原子类型
            solvent_atom_types = solvent_atom_types or [5, 6] # 默认溶剂原子类型
            
            self.logger.warning("使用默认原子类型，请检查是否合适")
        
        # 运行溶剂分布分析
        solvent_results = self.solvent_analyzer.analyze_solvent_distribution(
            dump_file, 
            atom_counts=atom_counts,
            cation_atom_types=cation_atom_types,
            anion_atom_types=anion_atom_types,
            solvent_atom_types=solvent_atom_types,
            cutoff_distance=cutoff_distance,
            cip_distance=cip_distance,
            ssip_distance=ssip_distance
        )
        
        if not solvent_results:
            self.logger.error("溶剂分布分析失败")
            return {'error': 'Solvent distribution analysis failed'}
        
        # 生成溶剂分布图
        bar_img = self.visualizer.plot_solvent_distribution(
            solvent_results, 
            f"solvent_distribution_{sample_name}_bar.png", 
            plot_type='bar'
        )
        
        pie_img = self.visualizer.plot_solvent_distribution(
            solvent_results, 
            f"solvent_distribution_{sample_name}_pie.png", 
            plot_type='pie'
        )
        
        heatmap_img = self.visualizer.plot_solvent_distribution(
            solvent_results, 
            f"solvent_distribution_{sample_name}_heatmap.png", 
            plot_type='heatmap'
        )
        
        # 生成离子对分析图表
        ion_pair_img = self.visualizer.plot_ion_pair_analysis(
            solvent_results, 
            f"ion_pair_analysis_{sample_name}.png"
        )
        
        # 将平均溶剂化数据导出为Excel
        avg_solvation_data = []
        for cation_id, solvents in solvent_results['avg_cation_solvation'].items():
            row = {'cation_id': cation_id}
            row.update(solvents)
            avg_solvation_data.append(row)
        
        # 获取离子对分析结果
        ion_pair_analysis = solvent_results.get('ion_pair_analysis', {})
        
        # 创建Excel文件
        excel_path = os.path.join(output_sample_dir, f'solvent_distribution_{sample_name}.xlsx')
        
        with pd.ExcelWriter(excel_path) as writer:
            # 溶剂化数据表
            if avg_solvation_data:
                df_solvation = pd.DataFrame(avg_solvation_data)
                df_solvation.to_excel(writer, sheet_name='溶剂化数据', index=False)
            
            # 离子对分析表
            ion_pair_data = {
                '类型': ['CIP (接触离子对)', 'SSIP (溶剂分离离子对)', 'AGG (聚集体)', '自由阳离子'],
                '数量': [
                    solvent_results.get('cip_count', 0),
                    solvent_results.get('ssip_count', 0),
                    solvent_results.get('agg_count', 0),
                    solvent_results.get('free_cation_count', 0)
                ],
                '百分比': [
                    solvent_results.get('cip_percentage', 0),
                    solvent_results.get('ssip_percentage', 0),
                    solvent_results.get('agg_percentage', 0),
                    solvent_results.get('free_cation_percentage', 0)
                ]
            }
            
            df_ion_pairs = pd.DataFrame(ion_pair_data)
            df_ion_pairs.to_excel(writer, sheet_name='离子对分析', index=False)
        
            # 溶剂化分布表
            solvation_dist = solvent_results.get('solvation_distribution', {})
            if solvation_dist:
                solvation_dist_data = {
                    '溶剂化数': list(solvation_dist.keys()),
                    '阳离子数': list(solvation_dist.values())
                }
                df_solvation_dist = pd.DataFrame(solvation_dist_data)
                df_solvation_dist = df_solvation_dist.sort_values('溶剂化数')
                df_solvation_dist.to_excel(writer, sheet_name='溶剂化分布', index=False)
        
        results = {
            'sample_name': sample_name,
            'dump_file': dump_file,
            'cutoff_distance': cutoff_distance,
            'cip_distance': cip_distance,
            'ssip_distance': ssip_distance,
            'cation_atom_types': cation_atom_types,
            'anion_atom_types': anion_atom_types,
            'solvent_atom_types': solvent_atom_types,
            'frames_analyzed': solvent_results.get('frames', 0),
            'cation_count': solvent_results.get('cation_count', 0),
            'anion_count': solvent_results.get('anion_count', 0),
            'solvent_count': solvent_results.get('solvent_count', 0),
            'cip_count': solvent_results.get('cip_count', 0),
            'ssip_count': solvent_results.get('ssip_count', 0),
            'agg_count': solvent_results.get('agg_count', 0),
            'free_cation_count': solvent_results.get('free_cation_count', 0),
            'bar_chart': bar_img,
            'pie_chart': pie_img,
            'heatmap': heatmap_img,
            'ion_pair_chart': ion_pair_img,
            'excel_file': excel_path
        }
        
        # 保存结果到JSON文件
        result_json_path = os.path.join(output_sample_dir, f'solvent_analysis_{sample_name}.json')
        with open(result_json_path, 'w') as f:
            # 过滤掉不可序列化的对象
            serializable_results = {k: v for k, v in results.items() 
                                   if isinstance(v, (str, int, float, list, dict, bool, type(None)))}
            json.dump(serializable_results, f, indent=2)
        
        self.logger.info(f"溶剂分析完成: {sample_name}")
        return results
    
    def run_comprehensive_analysis(self, sample_name: str) -> Dict[str, Any]:
        """
        运行综合分析，包括RDF、MSD和溶剂分析

        Args:
            sample_name: 样品名称

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"运行综合分析: {sample_name}")
        
        # 运行各个分析
        rdf_results = self.run_rdf_analysis(sample_name)
        msd_results = self.run_msd_analysis(sample_name)
        solvent_results = self.run_solvent_analysis(sample_name)
        
        # 创建综合可视化
        dashboard_data = {}
        
        if 'diffusion_coefficients' in msd_results:
            dashboard_data['diffusion'] = {sample_name: msd_results['diffusion_coefficients']}
        
        if 'rdf_image_path' in rdf_results:
            dashboard_data['rdf_data'] = rdf_results['rdf_image_path']
        
        if 'error' not in solvent_results:
            dashboard_data['solvation'] = {
                'avg_cation_solvation': solvent_results.get('avg_cation_solvation', {}),
                'avg_anion_solvation': solvent_results.get('avg_anion_solvation', {})
            }
        
        # 创建综合仪表板
        dashboard_img = self.visualizer.create_comparison_dashboard(
            dashboard_data,
            f"dashboard_{sample_name}.png"
        )
        
        results = {
            'sample_name': sample_name,
            'rdf_analysis': rdf_results,
            'msd_analysis': msd_results,
            'solvent_analysis': solvent_results,
            'dashboard': dashboard_img
        }
        
        # 保存结果到JSON文件
        output_sample_dir = os.path.join(self.output_dir, sample_name)
        os.makedirs(output_sample_dir, exist_ok=True)
        result_json_path = os.path.join(output_sample_dir, f'comprehensive_analysis_{sample_name}.json')
        
        with open(result_json_path, 'w') as f:
            # 过滤掉不可序列化的对象
            serializable_results = {
                'sample_name': sample_name,
                'dashboard': dashboard_img
            }
            json.dump(serializable_results, f, indent=2)
        
        self.logger.info(f"综合分析完成: {sample_name}")
        return results
    
    def run_batch_analysis(self, sample_names: List[str], analysis_types: List[str] = None) -> Dict[str, Any]:
        """
        批量运行分析任务

        Args:
            sample_names: 样品名称列表
            analysis_types: 分析类型列表，可选值: 'rdf', 'msd', 'solvent', 'comprehensive'，默认全部

        Returns:
            包含分析结果的字典
        """
        self.logger.info(f"批量运行分析任务: {sample_names}")
        
        # 默认运行所有类型的分析
        if analysis_types is None:
            analysis_types = ['rdf', 'msd', 'solvent', 'comprehensive', 'electrochemical']
        
        results = defaultdict(dict)
        
        for sample_name in sample_names:
            self.logger.info(f"处理样品: {sample_name}")
            
            for analysis_type in analysis_types:
                if analysis_type == 'rdf':
                    results['rdf'][sample_name] = self.run_rdf_analysis(sample_name)
                elif analysis_type == 'msd':
                    results['msd'][sample_name] = self.run_msd_analysis(sample_name)
                elif analysis_type == 'solvent':
                    results['solvent'][sample_name] = self.run_solvent_analysis(sample_name)
                elif analysis_type == 'comprehensive':
                    results['comprehensive'][sample_name] = self.run_comprehensive_analysis(sample_name)
                elif analysis_type == 'electrochemical':
                    # 对于电化学分析，我们需要额外的参数，这里使用默认值
                    results['electrochemical'][sample_name] = self.run_electrochemical_analysis(
                        sample_name=sample_name,
                        temperature=298.15,  # 室温 (K)
                        concentration=1.0    # 默认浓度 (mol/L)
                    )
                else:
                    self.logger.warning(f"未知的分析类型: {analysis_type}")
        
        # 如果有多个样品的MSD分析，添加比较分析
        if 'msd' in results and len(results['msd']) > 1:
            results['comparative_msd'] = self.run_comparative_msd_analysis(sample_names)
            
        # 如果有多个样品的电化学分析，添加比较分析
        if 'electrochemical' in results and len(results['electrochemical']) > 1:
            results['comparative_electrochemical'] = self.run_comparative_electrochemical_analysis(
                sample_names=sample_names,
                output_name="comparative_electrochemical",
                temperature=298.15,  # 室温 (K)
                concentrations=1.0   # 默认浓度 (mol/L)
            )
        
        # 保存批处理结果摘要
        summary = {
            'sample_names': sample_names,
            'analysis_types': analysis_types,
            'completion_status': {
                analysis_type: {
                    sample_name: 'error' not in results[analysis_type].get(sample_name, {})
                    for sample_name in sample_names
                } for analysis_type in analysis_types
            }
        }
        
        summary_path = os.path.join(self.output_dir, 'batch_analysis_summary.json')
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        self.logger.info(f"批量分析任务完成")
        return dict(results) 