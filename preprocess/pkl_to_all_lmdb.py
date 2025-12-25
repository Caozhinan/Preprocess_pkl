#!/usr/bin/env python3
"""将多个 PKL 文件合并为一个 LMDB，使用连续数字键，并在写入前检测 NaN（增强鲁棒性版）。"""

import argparse
import os
import pickle
from pathlib import Path

import lmdb
import numpy as np


def get_file_size(filepath: str) -> int:
    """Get file size in bytes."""
    return os.path.getsize(filepath) if os.path.exists(filepath) else 0


def has_nan_values(data: dict):
    """
    鲁棒地检测样本字典中是否包含 NaN。
    会对字典中的每个键值对进行深度递归检查。
    
    支持类型:
    - Numpy Arrays (检查 dtype 为 floating/complex 的)
    - PyTorch Tensors (尝试调用 .isnan() 或转 numpy)
    - Python Lists, Tuples, Sets (递归检查)
    - Nested Dictionaries (递归检查)
    - Standard Floats
    
    返回:
        has_nan: bool, 是否含 NaN
        bad_keys: set[str], 含有 NaN 的顶层键名集合
    """
    bad_keys = set()

    def check_recursive(value):
        """递归检查单个值是否包含 NaN"""
        try:
            # 1. 基础浮点数
            if isinstance(value, float):
                return np.isnan(value)

            # 2. Numpy 数组
            if isinstance(value, np.ndarray):
                # 只有浮点数或复数才可能有 NaN
                if np.issubdtype(value.dtype, np.floating) or np.issubdtype(value.dtype, np.complexfloating):
                    return np.isnan(value).any()
                # 如果是 Object 类型的数组（很少见但可能），递归检查元素
                elif value.dtype == object:
                    for item in value.flat:
                        if check_recursive(item):
                            return True
                return False

            # 3. 处理 PyTorch Tensor (Duck Typing，不强制依赖 torch 库)
            # 检查是否有 .isnan() 方法 (Tensor 特征)
            if hasattr(value, 'isnan'):
                try:
                    # item() 用于将单元素 tensor 转为 python bool，any() 用于多元素
                    res = value.isnan().any()
                    if hasattr(res, 'item'):
                        return res.item()
                    return bool(res)
                except:
                    pass
            
            # 兼容性回退：如果对象有 .numpy() 方法（如 Tensor），转为 numpy 再查
            if hasattr(value, 'numpy'):
                try:
                    arr = value.detach().cpu().numpy()
                    return check_recursive(arr)
                except:
                    pass

            # 4. 容器类型递归 (List, Tuple, Set)
            if isinstance(value, (list, tuple, set)):
                for item in value:
                    if check_recursive(item):
                        return True
                return False

            # 5. 字典递归
            if isinstance(value, dict):
                for v in value.values():
                    if check_recursive(v):
                        return True
                return False

        except Exception:
            # 如果在检测过程中发生任何无法预料的类型错误，保守起见视为无 NaN，或者你可以选择打印警告
            return False
        
        return False

    # 遍历顶层字典
    for key, value in data.items():
        if check_recursive(value):
            bad_keys.add(key)

    return (len(bad_keys) > 0), bad_keys


def main():
    parser = argparse.ArgumentParser(
        description="Merge PKL files to LMDB with numeric keys and robust NaN checking"
    )
    parser.add_argument("--csv_file", required=True, help="Input CSV file")
    parser.add_argument("--output_lmdb", required=True, help="Output LMDB path")
    args = parser.parse_args()

    # Collect PKL files and calculate total size
    pkl_files = []
    total_size = 0

    with open(args.csv_file, "r") as f:
        lines = f.readlines()

    # Skip header
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split(",")
        if len(parts) >= 3:
            receptor = parts[0].strip()
            name = parts[2].strip()

            # Construct PKL file path
            target_dir = os.path.dirname(receptor)
            pkl_path = os.path.join(target_dir, "output", f"{name}_features_with_masif.pkl")

            if os.path.exists(pkl_path):
                size = get_file_size(pkl_path)
                pkl_files.append((pkl_path, name, size))
                total_size += size
                print(f"Found: {pkl_path} ({size} bytes)")
            else:
                print(f"Missing: {pkl_path}")

    if not pkl_files:
        print("No valid PKL files found!")
        return

    # Calculate LMDB map size (1.2x total PKL size)
    map_size = int(total_size * 1.1)
    print(f"Total PKL size: {total_size / (1024**3):.2f} GB")
    print(f"LMDB map size: {map_size / (1024**3):.2f} GB")

    # Create LMDB database
    output_dir = os.path.dirname(args.output_lmdb)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    env = lmdb.open(
        args.output_lmdb,
        map_size=map_size,
        subdir=True,
        readonly=False,
        metasync=False,
        sync=False,
        map_async=True,
        writemap=True,
        meminit=False,
        max_readers=1,
    )

    count = 0
    total_data_size = 0
    skipped_files = []

    txn = env.begin(write=True)
    try:
        for pkl_path, name, size in pkl_files:
            try:
                with open(pkl_path, "rb") as f:
                    data = pickle.load(f)

                # 处理 list 格式
                if isinstance(data, list):
                    if len(data) > 0:
                        data = data[0]
                    else:
                        print(f"Error: Empty list in {pkl_path}")
                        skipped_files.append(name)
                        continue

                # 确保数据是字典格式
                if not isinstance(data, dict):
                    print(f"Error: Unexpected data type {type(data)} in {pkl_path}")
                    skipped_files.append(name)
                    continue

                # 检查 NaN：一旦任意键出现 NaN，整条样本跳过
                has_nan, bad_keys = has_nan_values(data)
                if has_nan:
                    bad_keys_str = ", ".join(sorted(bad_keys))
                    print(
                        f"Warning: NaN values detected in sample {name} "
                        f"(keys: {bad_keys_str}), skipping this sample."
                    )
                    skipped_files.append(name)
                    continue

                # 确保数据中包含原始名称
                if "pdbid" not in data:
                    data["pdbid"] = name

                # 使用连续数字键
                key = f"{count}".encode("ascii")
                value = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
                txn.put(key, value)

                total_data_size += len(value)
                count += 1

                # 每 1000 个样本提交一次事务
                if count % 1000 == 0:
                    txn.commit()
                    txn = env.begin(write=True)
                    print(f"  Committed at {count} samples")

            except Exception as e:
                print(f"Error processing {pkl_path}: {e}")
                skipped_files.append(name)

        # 提交最后一批数据
        txn.commit()

    except Exception as e:
        print(f"Fatal error during conversion: {e}")
        txn.abort()
        env.close()
        raise

    # 写入元数据 (在新的事务中)
    with env.begin(write=True) as meta_txn:
        metadata = {
            "total_samples": count,
            "skipped_files": skipped_files,
            "source_file": args.csv_file,
            "total_size_bytes": total_data_size,
        }
        meta_txn.put(b"__metadata__", pickle.dumps(metadata))

    env.close()

    print(f"\n✅ LMDB database created: {args.output_lmdb}")
    print(f"Total entries: {count}")
    print(f"Skipped files: {len(skipped_files)}")
    if skipped_files:
        print(f"Skipped (first 10): {skipped_files[:10]}")
    print(f"Data size: {total_data_size / (1024**3):.2f} GB")


if __name__ == "__main__":
    main()