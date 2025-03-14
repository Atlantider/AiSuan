#!/usr/bin/env python
"""
修复XYZ文件中的周期性边界问题的命令行工具
"""
import argparse
import sys
import os
from pathlib import Path

# 添加父目录到sys.path，使模块导入正常工作
sys.path.append(str(Path(__file__).resolve().parent.parent.parent))

from molyte_cursor.src.analysis.trajectory_processor import TrajectoryProcessor
from molyte_cursor.src.utils.logger import Logger

def main():
    """命令行入口函数"""
    # 配置参数解析器
    parser = argparse.ArgumentParser(
        description="修复XYZ文件中由于周期性边界条件导致的分子断裂问题"
    )
    
    # 添加参数
    parser.add_argument(
        "input",
        help="输入XYZ文件路径或包含XYZ文件的目录"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="输出XYZ文件路径（仅处理单个文件时有效）"
    )
    
    parser.add_argument(
        "-b", "--bond-threshold",
        type=float,
        default=2.0,
        help="判断原子间连接的距离阈值（埃），默认值：2.0"
    )
    
    parser.add_argument(
        "-w", "--wrap-only",
        action="store_true",
        help="仅将原子包装回盒子内，不进行分子修复"
    )
    
    parser.add_argument(
        "-p", "--pattern",
        default="*.xyz",
        help="处理目录时的文件匹配模式，默认为'*.xyz'"
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
    
    # 创建轨迹处理器
    processor = TrajectoryProcessor()
    
    input_path = Path(args.input)
    
    # 检查输入路径是否存在
    if not input_path.exists():
        logger.error(f"输入路径不存在: {input_path}")
        sys.exit(1)
    
    # 根据输入类型选择处理方式
    if input_path.is_file():
        # 处理单个文件
        output_path = args.output if args.output else None
        success = processor.fix_periodic_molecules(
            str(input_path),
            output_path,
            bond_threshold=args.bond_threshold,
            unwrap=not args.wrap_only
        )
        
        if success:
            logger.info("处理成功完成！")
        else:
            logger.error("处理失败！")
            sys.exit(1)
    
    elif input_path.is_dir():
        # 处理目录
        if args.output:
            logger.warning("处理目录时忽略输出文件参数，将在目录中创建fixed_trajectories子目录")
        
        count = processor.process_trajectory_directory(
            str(input_path),
            file_pattern=args.pattern,
            unwrap=not args.wrap_only
        )
        
        if count > 0:
            logger.info(f"成功处理了 {count} 个文件！")
        else:
            logger.error("没有成功处理任何文件！")
            sys.exit(1)
    
    else:
        logger.error(f"无法确定输入路径的类型: {input_path}")
        sys.exit(1)

if __name__ == "__main__":
    main() 