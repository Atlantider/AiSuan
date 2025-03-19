"""
电解液计算工作流

处理电解液计算的完整工作流程，包括输入解析、文件生成、计算提交和结果分析。
"""

import os
import sys
import time
import json
import logging
import subprocess
import tempfile
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

# 导入molyte_cursor模块
from molyte_cursor.src.io.inp_reader import read_inp_file, parse_inp_content, INPReader
from molyte_cursor.src.io.file_generator import LAMMPSFileGenerator, FileGenerator
from molyte_cursor.src.utils.command_executor import CommandExecutor
from molyte_cursor.src.analysis.analyzer import Analyzer
from molyte_cursor.src.utils.logger import setup_logger

class ElectrolyteWorkflow:
    """
    电解液计算工作流
    
    管理电解液计算的完整流程，包括：
    1. 解析INP输入文件
    2. 生成LAMMPS输入文件
    3. 提交LAMMPS计算
    4. 分析计算结果
    5. 生成输出报告
    """
    
    def __init__(self, 
                 inp_file_path: Optional[str] = None, 
                 inp_content: Optional[str] = None,
                 output_dir: Optional[str] = None,
                 log_level: str = 'INFO'):
        """
        初始化工作流
        
        Args:
            inp_file_path: INP文件路径（与inp_content二选一）
            inp_content: INP文件内容（与inp_file_path二选一）
            output_dir: 输出目录，默认为当前目录下的output/timestamp
            log_level: 日志级别
        """
        # 设置输出目录
        if not output_dir:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_dir = os.path.join(os.getcwd(), f"output_{timestamp}")
        
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 设置日志
        log_file = os.path.join(self.output_dir, 'workflow.log')
        self.logger = setup_logger(__name__, log_file, log_level)
        
        # 保存输入参数
        self.inp_file_path = inp_file_path
        self.inp_content = inp_content
        
        # 初始化状态变量
        self.config = None
        self.generated_files = []
        self.calculation_result = None
        self.status = 'initialized'
        
        # 创建工作目录
        self.work_dir = os.path.join(self.output_dir, 'work')
        os.makedirs(self.work_dir, exist_ok=True)
        
        # 初始化命令执行器
        self.cmd_executor = CommandExecutor(
            working_dir=self.work_dir,
            log_file=os.path.join(self.output_dir, 'commands.log')
        )
        
        self.logger.info(f"初始化电解液计算工作流, 输出目录: {self.output_dir}")
    
    def run(self) -> Dict[str, Any]:
        """
        运行完整工作流
        
        Returns:
            Dict: 计算结果字典
        """
        try:
            self.logger.info("开始执行电解液计算工作流")
            
            # 步骤1：解析INP输入
            self.parse_input()
            
            # 步骤2：生成输入文件
            self.generate_input_files()
            
            # 步骤3：运行LAMMPS计算
            self.run_calculation()
            
            # 步骤4：分析结果
            self.analyze_results()
            
            # 步骤5：生成报告
            self.generate_report()
            
            self.status = 'completed'
            self.logger.info("电解液计算工作流完成")
            
            return self.get_result_summary()
            
        except Exception as e:
            self.status = 'failed'
            self.logger.error(f"工作流执行失败: {str(e)}", exc_info=True)
            
            return {
                'status': 'failed',
                'error': str(e),
                'output_dir': self.output_dir
            }
    
    def parse_input(self) -> None:
        """解析INP输入，获取计算配置"""
        self.logger.info("步骤1: 解析INP输入")
        
        try:
            # 根据提供的输入类型选择解析方法
            if self.inp_file_path:
                self.logger.info(f"从文件读取INP: {self.inp_file_path}")
                self.config = read_inp_file(self.inp_file_path)
                
                # 复制输入文件到输出目录
                import shutil
                dest_path = os.path.join(self.output_dir, os.path.basename(self.inp_file_path))
                shutil.copy2(self.inp_file_path, dest_path)
                self.logger.info(f"INP文件已复制到: {dest_path}")
                
            elif self.inp_content:
                self.logger.info("从字符串解析INP内容")
                self.config = parse_inp_content(self.inp_content)
                
                # 保存内容到输入文件
                inp_path = os.path.join(self.output_dir, 'input.inp')
                with open(inp_path, 'w', encoding='utf-8') as f:
                    f.write(self.inp_content)
                self.logger.info(f"INP内容已保存到: {inp_path}")
                
            else:
                raise ValueError("未提供INP文件路径或内容")
            
            # 保存解析后的配置
            config_path = os.path.join(self.output_dir, 'config.json')
            with open(config_path, 'w', encoding='utf-8') as f:
                json.dump(self.config, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"配置已解析并保存到: {config_path}")
            self.logger.debug(f"解析配置: {json.dumps(self.config, indent=2)}")
            
            # 更新状态
            self.status = 'input_parsed'
            
        except Exception as e:
            self.logger.error(f"解析INP输入失败: {str(e)}", exc_info=True)
            raise ValueError(f"解析INP输入失败: {str(e)}")
    
    def generate_input_files(self) -> None:
        """生成LAMMPS输入文件"""
        self.logger.info("步骤2: 生成LAMMPS输入文件")
        
        try:
            # 确保配置已解析
            if not self.config:
                raise ValueError("配置未解析，请先调用parse_input()方法")
            
            # 使用FileGenerator生成文件
            file_generator = FileGenerator(
                working_dir=self.work_dir,
                temp_dir=os.path.join(self.output_dir, 'temp')
            )
            
            # 生成文件
            self.generated_files = file_generator.generate_files(self.config)
            
            # 记录生成的文件
            files_list_path = os.path.join(self.output_dir, 'generated_files.json')
            with open(files_list_path, 'w', encoding='utf-8') as f:
                json.dump(self.generated_files, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"已生成{len(self.generated_files)}个输入文件")
            self.logger.debug(f"生成的文件: {json.dumps(self.generated_files, indent=2)}")
            
            # 更新状态
            self.status = 'files_generated'
            
        except Exception as e:
            self.logger.error(f"生成输入文件失败: {str(e)}", exc_info=True)
            raise ValueError(f"生成输入文件失败: {str(e)}")
    
    def run_calculation(self) -> None:
        """运行LAMMPS计算"""
        self.logger.info("步骤3: 运行LAMMPS计算")
        
        try:
            # 确保文件已生成
            if not self.generated_files:
                raise ValueError("输入文件未生成，请先调用generate_input_files()方法")
            
            # 找到主要的LAMMPS输入文件
            main_input_file = None
            for file_info in self.generated_files:
                if file_info.get('type') == 'lammps_input' and file_info.get('is_main', False):
                    main_input_file = file_info.get('path')
                    break
            
            if not main_input_file:
                # 尝试查找任何LAMMPS输入文件
                for file_info in self.generated_files:
                    if file_info.get('type') == 'lammps_input':
                        main_input_file = file_info.get('path')
                        break
            
            if not main_input_file:
                raise ValueError("找不到LAMMPS输入文件")
            
            self.logger.info(f"使用LAMMPS输入文件: {main_input_file}")
            
            # 执行LAMMPS计算
            start_time = time.time()
            cmd = f"lammps -in {main_input_file}"
            
            self.logger.info(f"执行命令: {cmd}")
            result = self.cmd_executor.execute(cmd)
            
            end_time = time.time()
            duration = end_time - start_time
            
            # 保存命令输出
            output_path = os.path.join(self.output_dir, 'lammps_output.log')
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(result.stdout)
            
            self.logger.info(f"LAMMPS计算完成，耗时: {duration:.2f}秒")
            self.logger.info(f"输出已保存到: {output_path}")
            
            # 检查计算是否成功
            if result.returncode != 0:
                error_path = os.path.join(self.output_dir, 'lammps_error.log')
                with open(error_path, 'w', encoding='utf-8') as f:
                    f.write(result.stderr)
                
                self.logger.error(f"LAMMPS计算失败，错误码: {result.returncode}")
                self.logger.error(f"错误信息已保存到: {error_path}")
                raise RuntimeError(f"LAMMPS计算失败，错误码: {result.returncode}")
            
            # 更新状态
            self.status = 'calculation_completed'
            
        except Exception as e:
            self.logger.error(f"运行LAMMPS计算失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"运行LAMMPS计算失败: {str(e)}")
    
    def analyze_results(self) -> None:
        """分析计算结果"""
        self.logger.info("步骤4: 分析计算结果")
        
        try:
            # 确保计算已完成
            if self.status != 'calculation_completed':
                raise ValueError("计算未完成，请先调用run_calculation()方法")
            
            # 创建分析器
            analyzer = Analyzer(
                base_dir=self.work_dir,
                output_dir=os.path.join(self.output_dir, 'analysis')
            )
            
            # 提取计算类型
            calculation_types = self.config.get('calculation_types', ['conductivity'])
            
            # 运行相应的分析
            analysis_results = {}
            
            # RDF分析
            if 'rdf' in calculation_types or 'conductivity' in calculation_types:
                self.logger.info("执行RDF分析")
                rdf_results = analyzer.run_rdf_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['rdf'] = rdf_results
            
            # MSD分析
            if 'msd' in calculation_types or 'conductivity' in calculation_types:
                self.logger.info("执行MSD分析")
                msd_results = analyzer.run_msd_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['msd'] = msd_results
            
            # 电导率分析
            if 'conductivity' in calculation_types:
                self.logger.info("执行电导率分析")
                conductivity_results = analyzer.run_electrochemical_analysis(
                    sample_name=self.config.get('formulation_name', 'sample'),
                    temperature=self.config.get('temperature', 298),
                    concentration=self.config.get('concentration', 1.0)
                )
                analysis_results['conductivity'] = conductivity_results
            
            # 溶剂分析
            if 'solvent' in calculation_types:
                self.logger.info("执行溶剂分析")
                solvent_results = analyzer.run_solvent_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['solvent'] = solvent_results
            
            # 综合分析
            if 'comprehensive' in calculation_types:
                self.logger.info("执行综合分析")
                comprehensive_results = analyzer.run_comprehensive_analysis(
                    sample_name=self.config.get('formulation_name', 'sample')
                )
                analysis_results['comprehensive'] = comprehensive_results
            
            # 保存分析结果
            self.calculation_result = analysis_results
            
            # 将结果保存到文件
            results_path = os.path.join(self.output_dir, 'analysis_results.json')
            with open(results_path, 'w', encoding='utf-8') as f:
                # 将复杂对象转换为简单字典
                serializable_results = self._make_serializable(analysis_results)
                json.dump(serializable_results, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"分析结果已保存到: {results_path}")
            
            # 更新状态
            self.status = 'analysis_completed'
            
        except Exception as e:
            self.logger.error(f"分析计算结果失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"分析计算结果失败: {str(e)}")
    
    def generate_report(self) -> None:
        """生成结果报告"""
        self.logger.info("步骤5: 生成结果报告")
        
        try:
            # 确保分析已完成
            if self.status != 'analysis_completed':
                raise ValueError("分析未完成，请先调用analyze_results()方法")
            
            # 生成摘要报告
            summary = self.get_result_summary()
            
            # 保存报告
            report_path = os.path.join(self.output_dir, 'report.json')
            with open(report_path, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"结果报告已保存到: {report_path}")
            
            # 更新状态
            self.status = 'report_generated'
            
        except Exception as e:
            self.logger.error(f"生成结果报告失败: {str(e)}", exc_info=True)
            raise RuntimeError(f"生成结果报告失败: {str(e)}")
    
    def get_result_summary(self) -> Dict[str, Any]:
        """
        获取计算结果摘要
        
        Returns:
            Dict: 结果摘要字典
        """
        summary = {
            'status': self.status,
            'formulation_name': self.config.get('formulation_name', 'Unnamed') if self.config else 'Unknown',
            'output_dir': self.output_dir,
            'timestamp': datetime.now().isoformat(),
        }
        
        # 如果有计算结果，提取关键信息
        if self.calculation_result:
            # 从conductivity结果中提取关键数据
            if 'conductivity' in self.calculation_result:
                cond_result = self.calculation_result['conductivity']
                if isinstance(cond_result, dict):
                    summary['conductivity_mS_cm'] = cond_result.get('conductivity_mS_cm')
                    summary['molar_conductivity'] = cond_result.get('molar_conductivity')
                    summary['transference_number'] = cond_result.get('transference_number')
                    summary['viscosity_mPa_s'] = cond_result.get('viscosity_mPa_s')
                    
                    # 提取扩散系数
                    if 'diffusion_coeffs' in cond_result:
                        summary['diffusion_coefficients'] = cond_result['diffusion_coeffs']
            
            # 从comprehensive结果中提取其他信息
            if 'comprehensive' in self.calculation_result:
                comp_result = self.calculation_result['comprehensive']
                if isinstance(comp_result, dict) and 'density' not in summary:
                    summary['density'] = comp_result.get('density')
        
        return summary
    
    def _make_serializable(self, obj: Any) -> Any:
        """
        将对象转换为可序列化的格式
        
        Args:
            obj: 要转换的对象
            
        Returns:
            可序列化的对象
        """
        if isinstance(obj, dict):
            return {k: self._make_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_serializable(item) for item in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        else:
            # 尝试转换为字典
            try:
                return self._make_serializable(obj.__dict__)
            except AttributeError:
                # 如果无法转换，则转为字符串
                return str(obj)


def run_from_inp_file(inp_file_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    从INP文件运行完整电解液计算工作流
    
    Args:
        inp_file_path: INP文件路径
        output_dir: 输出目录
        
    Returns:
        Dict: 计算结果摘要
    """
    workflow = ElectrolyteWorkflow(inp_file_path=inp_file_path, output_dir=output_dir)
    return workflow.run()

def run_from_inp_content(inp_content: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    从INP内容字符串运行完整电解液计算工作流
    
    Args:
        inp_content: INP文件内容
        output_dir: 输出目录
        
    Returns:
        Dict: 计算结果摘要
    """
    workflow = ElectrolyteWorkflow(inp_content=inp_content, output_dir=output_dir)
    return workflow.run() 