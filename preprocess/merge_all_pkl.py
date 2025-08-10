#!/usr/bin/env python3
"""
merge_all_pkl.py - 合并多个复合物的PKL文件（支持多核并行）

用法:
python merge_all_pkl.py --input_dir /path/to/directory --output_file merged_complexes.pkl
python merge_all_pkl.py --input_pattern "/path/to/*_features_with_masif.pkl" --output_file merged.pkl
python merge_all_pkl.py --input_file_list my_pkls.txt --output_file merged.pkl
"""

import argparse
import pickle
import glob
from pathlib import Path
from typing import List, Dict, Any
import sys
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_pkl_file(pkl_path: str) -> List[Dict[Any, Any]]:
    """
    加载单个PKL文件

    Args:
        pkl_path: PKL文件路径

    Returns:
        复合物数据列表
    """
    try:
        with open(pkl_path, 'rb') as f:
            data = pickle.load(f)

        # 确保数据是列表格式
        if isinstance(data, list):
            return data
        else:
            # 如果是单个复合物，包装成列表
            return [data]

    except Exception as e:
        print(f"Error loading {pkl_path}: {e}")
        return []

def merge_pkl_files(pkl_files: List[str], output_file: str, num_workers: int = 12) -> None:
    """
    并行合并多个PKL文件到单个文件

    Args:
        pkl_files: PKL文件路径列表
        output_file: 输出文件路径
        num_workers: 并行核数
    """
    all_complexes = []

    print(f"Found {len(pkl_files)} PKL files to merge (using {num_workers} workers):")

    # 并发读取
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_file = {executor.submit(load_pkl_file, pkl_file): pkl_file for pkl_file in pkl_files}
        for i, future in enumerate(as_completed(future_to_file), 1):
            pkl_file = future_to_file[future]
            try:
                complexes = future.result()
                if complexes:
                    all_complexes.extend(complexes)
                    print(f"[{i}/{len(pkl_files)}] Added {len(complexes)} complex(es) from {pkl_file}")
                else:
                    print(f"[{i}/{len(pkl_files)}] Warning: No data loaded from {pkl_file}")
            except Exception as e:
                print(f"[{i}/{len(pkl_files)}] Error reading {pkl_file}: {e}")

    if not all_complexes:
        print("Error: No complexes found to merge!")
        sys.exit(1)

    print(f"\nTotal complexes to save: {len(all_complexes)}")

    # 保存合并后的数据
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(all_complexes, f)
        print(f"Successfully saved merged data to: {output_file}")

        # 验证保存的文件
        with open(output_file, 'rb') as f:
            verification_data = pickle.load(f)
        print(f"Verification: Loaded {len(verification_data)} complexes from saved file")

    except Exception as e:
        print(f"Error saving merged file: {e}")
        sys.exit(1)

def find_pkl_files(input_dir: str = None, input_pattern: str = None, input_file_list: str = None) -> List[str]:
    """
    查找PKL文件

    Args:
        input_dir: 输入目录路径
        input_pattern: 文件模式匹配
        input_file_list: 文件列表

    Returns:
        PKL文件路径列表
    """
    pkl_files = []

    if input_file_list:
        if not os.path.isfile(input_file_list):
            print(f"Error: input_file_list file not found: {input_file_list}")
            sys.exit(1)
        with open(input_file_list) as fin:
            pkl_files = [line.strip() for line in fin if line.strip()]
    elif input_pattern:
        pkl_files = glob.glob(input_pattern, recursive=True)
    elif input_dir:
        pattern = str(Path(input_dir) / "**" / "*_features_with_masif.pkl")
        pkl_files = glob.glob(pattern, recursive=True)
        if not pkl_files:
            pattern = str(Path(input_dir) / "**" / "*_features.pkl")
            pkl_files = glob.glob(pattern, recursive=True)

    # 只保留存在的文件
    pkl_files = [f for f in pkl_files if os.path.isfile(f)]
    return sorted(pkl_files)

def main():
    parser = argparse.ArgumentParser(
        description="合并多个复合物的PKL文件到单个文件（支持多核并行）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 合并目录下所有PKL文件
  python merge_all_pkl.py --input_dir /path/to/data --output_file merged.pkl

  # 使用通配符模式
  python merge_all_pkl.py --input_pattern "/path/to/*_features_with_masif.pkl" --output_file merged.pkl

  # 使用文件列表
  python merge_all_pkl.py --input_file_list my_pkls.txt --output_file merged.pkl

  # 合并当前目录下的所有PKL文件
  python merge_all_pkl.py --input_dir . --output_file all_complexes.pkl

  # 指定并行核数
  python merge_all_pkl.py --input_dir . --output_file all_complexes.pkl --num_workers 16
        """
    )

    # 输入选项（三选一）
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input_dir',
        type=str,
        help='输入目录路径，将递归查找所有 *_features_with_masif.pkl 文件'
    )
    input_group.add_argument(
        '--input_pattern',
        type=str,
        help='文件匹配模式，如 "/path/to/*_features_with_masif.pkl"'
    )
    input_group.add_argument(
        '--input_file_list',
        type=str,
        help='包含所有要合并的pkl文件路径的文本文件，每行一个路径'
    )

    # 输出文件
    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help='输出PKL文件路径'
    )

    # 可选参数
    parser.add_argument(
        '--dry_run',
        action='store_true',
        help='仅显示将要处理的文件，不实际合并'
    )

    parser.add_argument(
        '--num_workers',
        type=int,
        default=12,
        help='并行核数（默认12）'
    )

    args = parser.parse_args()

    # 查找PKL文件
    pkl_files = find_pkl_files(
        input_dir=getattr(args, 'input_dir', None),
        input_pattern=getattr(args, 'input_pattern', None),
        input_file_list=getattr(args, 'input_file_list', None)
    )

    if not pkl_files:
        print("Error: No PKL files found!")
        if args.input_dir:
            print(f"Searched in directory: {args.input_dir}")
        if args.input_pattern:
            print(f"Searched with pattern: {args.input_pattern}")
        if args.input_file_list:
            print(f"Checked file list: {args.input_file_list}")
        sys.exit(1)

    if args.dry_run:
        print(f"Found {len(pkl_files)} PKL files:")
        for pkl_file in pkl_files:
            print(f"  {pkl_file}")
        print(f"Would save merged data to: {args.output_file}")
        return

    # 执行合并
    merge_pkl_files(pkl_files, args.output_file, num_workers=args.num_workers)

if __name__ == "__main__":
    main()