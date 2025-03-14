#!/usr/bin/env python
"""
Molyte集成工作流模块，提供预处理和后处理的统一入口
"""
import os
import sys
import argparse
import shutil
from pathlib import Path
import subprocess
import datetime

from ..utils.logger import Logger
from ..io.excel_reader import ExcelReader
from ..simulation.simulator import Simulator
from ..analysis.gaussian_processor import GaussianProcessor
from ..analysis.trajectory_processor import TrajectoryProcessor

class IntegratedWorkflow:
    """Molyte集成工作流类"""
    
    def __init__(self, workspace_root=None):
        """初始化集成工作流
        
        Args:
            workspace_root: 工作空间根目录，默认为当前目录
        """
        self.logger = Logger().get_logger()
        
        # 设置工作空间根目录
        if workspace_root:
            self.workspace_root = Path(workspace_root)
        else:
            self.workspace_root = Path.cwd()
        
        # 创建标准工作目录结构
        self.lammps_workspace = self.workspace_root / "Lammps_Workspace"
        self.molecule_workspace = self.workspace_root / "Molecule_Workspace"
        self.lammps_results = self.workspace_root / "Lammps_Results"
        self.gaussian_results = self.workspace_root / "Gaussian_Results"
        
        # 创建所有工作目录
        self.lammps_workspace.mkdir(parents=True, exist_ok=True)
        self.molecule_workspace.mkdir(parents=True, exist_ok=True)
        self.lammps_results.mkdir(parents=True, exist_ok=True)
        self.gaussian_results.mkdir(parents=True, exist_ok=True)
        
        # 初始化组件
        self.excel_reader = ExcelReader()
        self.simulator = Simulator()
        self.gaussian_processor = GaussianProcessor()
        self.trajectory_processor = TrajectoryProcessor()
        
        # 记录已处理的系统
        self.processed_systems = []
    
    def preprocess(self, excel_file, run_lammps=True, run_gaussian=False, 
                  gaussian_method="B3LYP", gaussian_basis="6-31G(d)", 
                  gaussian_calc_type="opt freq"):
        """预处理阶段：设置并可选运行LAMMPS和分子计算
        
        Args:
            excel_file: Excel输入文件路径
            run_lammps: 是否运行LAMMPS模拟
            run_gaussian: 是否运行高斯计算
            gaussian_method: 高斯计算方法
            gaussian_basis: 高斯计算基组
            gaussian_calc_type: 高斯计算类型
            
        Returns:
            处理的系统列表
        """
        self.logger.info(f"开始预处理阶段，输入文件：{excel_file}")
        
        # 读取Excel文件
        systems = self.excel_reader.read_excel(excel_file)
        if not systems:
            self.logger.error(f"无法从Excel文件 {excel_file} 读取系统")
            return []
        
        # 设置LAMMPS工作目录
        # 复制Excel文件到LAMMPS工作空间，添加时间戳
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        excel_basename = Path(excel_file).name
        lammps_excel_copy = self.lammps_workspace / f"{timestamp}_{excel_basename}"
        shutil.copy2(excel_file, lammps_excel_copy)
        
        # 处理每个系统
        for system in systems:
            self.logger.info(f"处理系统：{system.name}")
            
            # 记录系统信息
            system_info = {
                "name": system.name,
                "timestamp": timestamp,
                "temperature": getattr(system, 'temperature', 300.0),
                "excel_file": str(lammps_excel_copy),
                "lammps_dir": None,
                "gaussian_dir": None
            }
            
            # 设置并运行LAMMPS模拟
            if run_lammps:
                self.logger.info(f"设置LAMMPS模拟：{system.name}")
                
                # 设置模拟，如果submit=True则会提交作业
                simulation_success = self.simulator.simulate(system, submit_job=True)
                
                if simulation_success:
                    # 记录LAMMPS目录
                    lammps_dir = self._get_lammps_dir(system)
                    system_info["lammps_dir"] = str(lammps_dir)
                    self.logger.info(f"LAMMPS模拟目录：{lammps_dir}")
                else:
                    self.logger.warning(f"LAMMPS模拟设置失败：{system.name}")
            
            # 设置并运行分子计算
            if run_gaussian:
                self.logger.info(f"设置分子计算：{system.name}")
                
                # 为该系统创建高斯工作目录
                gaussian_dir = self.molecule_workspace / f"{timestamp}_{system.name}"
                gaussian_dir.mkdir(parents=True, exist_ok=True)
                
                # 将系统信息写入单独的Excel文件
                system_excel = gaussian_dir / f"{system.name}.xlsx"
                # TODO: 实现将单个系统写入Excel文件的功能
                # 这里简单复制原Excel文件
                shutil.copy2(excel_file, system_excel)
                
                # 生成高斯输入文件
                summary_df, plot_files = self.gaussian_processor.calculate_and_visualize_properties(
                    system_excel,
                    gaussian_dir,
                    run_gaussian=True,  # 运行高斯
                    method=gaussian_method,
                    basis_set=gaussian_basis,
                    calc_type=gaussian_calc_type
                )
                
                # 记录高斯目录
                system_info["gaussian_dir"] = str(gaussian_dir)
            
            # 添加到已处理系统列表
            self.processed_systems.append(system_info)
        
        self.logger.info(f"预处理阶段完成，处理了 {len(systems)} 个系统")
        return self.processed_systems
    
    def postprocess(self, run_lammps_analysis=True, run_gaussian_analysis=True):
        """后处理阶段：处理LAMMPS和分子计算结果
        
        Args:
            run_lammps_analysis: 是否运行LAMMPS结果分析
            run_gaussian_analysis: 是否运行分子计算结果分析
            
        Returns:
            后处理结果字典
        """
        self.logger.info("开始后处理阶段")
        
        results = {
            "lammps": [],
            "gaussian": []
        }
        
        # 获取所有LAMMPS和Gaussian目录
        if not self.processed_systems:
            # 如果没有预处理过的系统，检查工作目录
            self.logger.info("没有预处理系统信息，扫描工作目录")
            
            # 扫描LAMMPS工作空间
            lammps_dirs = [d for d in self.lammps_workspace.iterdir() if d.is_dir()]
            
            # 扫描分子计算工作空间
            gaussian_dirs = [d for d in self.molecule_workspace.iterdir() if d.is_dir()]
            
            # 创建系统信息
            for lammps_dir in lammps_dirs:
                system_name = lammps_dir.name
                self.processed_systems.append({
                    "name": system_name,
                    "lammps_dir": str(lammps_dir),
                    "gaussian_dir": None
                })
            
            # 匹配高斯目录
            for gaussian_dir in gaussian_dirs:
                # 尝试从目录名提取系统名
                for system in self.processed_systems:
                    if system["name"] in gaussian_dir.name:
                        system["gaussian_dir"] = str(gaussian_dir)
                        break
                else:
                    # 如果没有匹配到已知系统，添加新系统
                    self.processed_systems.append({
                        "name": gaussian_dir.name,
                        "lammps_dir": None,
                        "gaussian_dir": str(gaussian_dir)
                    })
        
        # 处理每个系统
        for system_info in self.processed_systems:
            system_name = system_info["name"]
            self.logger.info(f"后处理系统：{system_name}")
            
            # 处理LAMMPS结果
            if run_lammps_analysis and system_info["lammps_dir"]:
                lammps_dir = Path(system_info["lammps_dir"])
                
                if lammps_dir.exists():
                    self.logger.info(f"处理LAMMPS结果：{lammps_dir}")
                    
                    # 创建系统结果目录
                    system_results_dir = self.lammps_results / system_name
                    system_results_dir.mkdir(parents=True, exist_ok=True)
                    
                    # 处理XYZ轨迹文件
                    xyz_files = list(lammps_dir.glob("*.xyz"))
                    for xyz_file in xyz_files:
                        self.logger.info(f"处理轨迹文件：{xyz_file}")
                        
                        # 修复周期性边界问题
                        fixed_xyz = system_results_dir / f"{xyz_file.stem}_fixed.xyz"
                        success = self.trajectory_processor.fix_periodic_molecules(
                            str(xyz_file),
                            str(fixed_xyz)
                        )
                        
                        if success:
                            results["lammps"].append({
                                "system": system_name,
                                "original_file": str(xyz_file),
                                "processed_file": str(fixed_xyz)
                            })
                    
                    # 复制其他关键结果文件
                    for pattern in ["*.log", "*.dat", "*.out"]:
                        for result_file in lammps_dir.glob(pattern):
                            dest_file = system_results_dir / result_file.name
                            if not dest_file.exists():
                                shutil.copy2(result_file, dest_file)
                                self.logger.debug(f"复制结果文件：{result_file} -> {dest_file}")
                
                else:
                    self.logger.warning(f"LAMMPS目录不存在：{lammps_dir}")
            
            # 处理高斯计算结果
            if run_gaussian_analysis and system_info["gaussian_dir"]:
                gaussian_dir = Path(system_info["gaussian_dir"])
                
                if gaussian_dir.exists():
                    self.logger.info(f"处理分子计算结果：{gaussian_dir}")
                    
                    # 创建系统结果目录
                    system_results_dir = self.gaussian_results / system_name
                    system_results_dir.mkdir(parents=True, exist_ok=True)
                    
                    # 解析高斯输出文件
                    gaussian_results, summary_df = self.gaussian_processor.parse_all_outputs(gaussian_dir)
                    
                    if gaussian_results:
                        # 生成可视化图表
                        plot_files = self.gaussian_processor.visualize_results(
                            gaussian_results,
                            system_results_dir
                        )
                        
                        # 复制结果摘要
                        if summary_df is not None:
                            summary_file = system_results_dir / "gaussian_results_summary.csv"
                            summary_df.to_csv(summary_file, index=False)
                            self.logger.info(f"保存结果摘要：{summary_file}")
                        
                        # 复制其他结果文件
                        for log_file in gaussian_dir.glob("*.log"):
                            dest_file = system_results_dir / log_file.name
                            if not dest_file.exists():
                                shutil.copy2(log_file, dest_file)
                        
                        results["gaussian"].append({
                            "system": system_name,
                            "results_dir": str(system_results_dir),
                            "plots": [str(p) for p in plot_files]
                        })
                
                else:
                    self.logger.warning(f"高斯目录不存在：{gaussian_dir}")
        
        self.logger.info(f"后处理阶段完成，处理了 {len(results['lammps'])} 个LAMMPS结果和 {len(results['gaussian'])} 个高斯结果")
        return results
    
    def _get_lammps_dir(self, system):
        """获取LAMMPS目录路径
        
        Args:
            system: 系统对象
            
        Returns:
            LAMMPS目录路径
        """
        # 从模拟器的配置中获取文件生成路径
        file_generate_path = self.simulator.file_generate_path
        
        # 生成系统目录名（包含温度）
        temperature = getattr(system, 'temperature', 300.0)
        system_dir_name = f"{system.name}_{int(temperature)}K"
        
        return file_generate_path / system_dir_name


def main():
    """命令行入口函数"""
    # 配置参数解析器
    parser = argparse.ArgumentParser(
        description="Molyte集成工作流：统一预处理和后处理流程"
    )
    
    # 添加模式选择参数
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--preprocess",
        action="store_true",
        help="执行预处理阶段（设置并运行LAMMPS和分子计算）"
    )
    mode_group.add_argument(
        "--postprocess",
        action="store_true",
        help="执行后处理阶段（处理LAMMPS和分子计算结果）"
    )
    
    # 添加通用参数
    parser.add_argument(
        "-w", "--workspace",
        help="工作空间根目录，默认为当前目录"
    )
    
    # 添加预处理参数
    parser.add_argument(
        "-i", "--input",
        help="输入Excel文件路径（预处理阶段必需）"
    )
    
    parser.add_argument(
        "--no-lammps",
        action="store_true",
        help="不运行LAMMPS模拟"
    )
    
    parser.add_argument(
        "--run-gaussian",
        action="store_true",
        help="运行高斯计算"
    )
    
    parser.add_argument(
        "--gaussian-method",
        default="B3LYP",
        help="高斯计算方法，默认为'B3LYP'"
    )
    
    parser.add_argument(
        "--gaussian-basis",
        default="6-31G(d)",
        help="高斯计算基组，默认为'6-31G(d)'"
    )
    
    parser.add_argument(
        "--gaussian-calc-type",
        default="opt freq",
        help="高斯计算类型，默认为'opt freq'"
    )
    
    # 添加后处理参数
    parser.add_argument(
        "--no-lammps-analysis",
        action="store_true",
        help="不运行LAMMPS结果分析"
    )
    
    parser.add_argument(
        "--no-gaussian-analysis",
        action="store_true",
        help="不运行分子计算结果分析"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="输出详细日志信息"
    )
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置日志级别
    logger = Logger(log_level="DEBUG" if args.verbose else "INFO").get_logger()
    
    # 创建集成工作流
    workflow = IntegratedWorkflow(workspace_root=args.workspace)
    
    # 根据模式执行不同操作
    if args.preprocess:
        # 检查输入文件是否提供
        if not args.input:
            logger.error("预处理阶段需要提供输入Excel文件")
            parser.print_help()
            sys.exit(1)
        
        # 执行预处理阶段
        processed_systems = workflow.preprocess(
            args.input,
            run_lammps=not args.no_lammps,
            run_gaussian=args.run_gaussian,
            gaussian_method=args.gaussian_method,
            gaussian_basis=args.gaussian_basis,
            gaussian_calc_type=args.gaussian_calc_type
        )
        
        if processed_systems:
            logger.info(f"预处理完成，处理了 {len(processed_systems)} 个系统")
            for system in processed_systems:
                logger.info(f"系统：{system['name']}")
                if system['lammps_dir']:
                    logger.info(f"  LAMMPS目录：{system['lammps_dir']}")
                if system['gaussian_dir']:
                    logger.info(f"  高斯目录：{system['gaussian_dir']}")
        else:
            logger.warning("预处理阶段未处理任何系统")
    
    elif args.postprocess:
        # 执行后处理阶段
        results = workflow.postprocess(
            run_lammps_analysis=not args.no_lammps_analysis,
            run_gaussian_analysis=not args.no_gaussian_analysis
        )
        
        if results["lammps"] or results["gaussian"]:
            logger.info("后处理完成")
            if results["lammps"]:
                logger.info(f"处理了 {len(results['lammps'])} 个LAMMPS结果")
                logger.info(f"结果保存在：{workflow.lammps_results}")
            if results["gaussian"]:
                logger.info(f"处理了 {len(results['gaussian'])} 个高斯结果")
                logger.info(f"结果保存在：{workflow.gaussian_results}")
        else:
            logger.warning("后处理阶段未处理任何结果")


if __name__ == "__main__":
    main() 