"""
主程序入口模块，协调整个程序的运行流程
"""
import sys
import argparse
from pathlib import Path
from ..utils.logger import Logger
from ..utils.config_loader import ConfigLoader
from ..io.excel_reader import ExcelReader
from ..simulation.simulator import Simulator

def parse_arguments():
    """解析命令行参数
    
    Returns:
        解析后的参数对象
    """
    parser = argparse.ArgumentParser(description="Molyte - 分子动力学模拟工具")
    
    parser.add_argument(
        "--input", "-i",
        type=str,
        help="输入Excel文件路径"
    )
    
    parser.add_argument(
        "--submit", "-s",
        action="store_true",
        help="是否提交计算作业"
    )
    
    parser.add_argument(
        "--log", "-l",
        type=str,
        help="日志文件路径"
    )
    
    parser.add_argument(
        "--config", "-c",
        type=str,
        help="配置目录路径"
    )
    
    return parser.parse_args()

def main():
    """主程序入口函数"""
    # 解析命令行参数
    args = parse_arguments()
    
    # 初始化日志
    logger = Logger(log_file=args.log)
    logger.info("Starting Molyte")
    
    try:
        # 加载配置
        config = ConfigLoader(args.config)
        config.load_paths_config()
        config.load_simulation_params()
        
        # 获取输入文件路径
        input_file = args.input
        if not input_file:
            input_file = config.get_file_path("input_xlsx")
            if not input_file:
                logger.error("No input file specified and no default input file found in config.")
                sys.exit(1)
        
        # 获取现有分子文件路径
        existing_molecules_file = config.get_file_path("existing_molecules")
        
        # 初始化Excel读取器
        excel_reader = ExcelReader(existing_molecules_file)
        
        # 读取Excel文件
        logger.info(f"Reading input file: {input_file}")
        systems = excel_reader.read_excel(input_file)
        logger.info(f"Found {len(systems)} systems in the input file.")
        
        # 初始化模拟器
        simulator = Simulator()
        
        # 处理每个系统
        successful_count = 0
        for i, system in enumerate(systems, start=1):
            logger.info(f"Processing system {i}/{len(systems)}: {system.name}")
            
            # 执行模拟
            success = simulator.simulate(system, args.submit)
            
            if success:
                successful_count += 1
                logger.info(f"Successfully processed system: {system.name}")
            else:
                logger.warning(f"Failed to process system: {system.name}")
        
        # 输出统计信息
        logger.info(f"Processed {len(systems)} systems, {successful_count} succeeded, {len(systems) - successful_count} failed.")
        logger.info("Molyte completed.")
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main() 