#!/usr/bin/env python3  
"""  
清理临时文件并合并所有PKL文件  
"""  
  
import argparse  
import pickle  
import os  
import sys  
import shutil  
import glob  
from pathlib import Path  
from concurrent.futures import ProcessPoolExecutor  
import traceback  
  
def cleanup_directory(target_dir, name):  
    """  
    清理单个目录，只保留必要的PKL文件  
    """  
    try:  
        output_dir = os.path.join(target_dir, "output")  
        if not os.path.exists(output_dir):  
            return False, f"Output directory not found: {output_dir}"  
              
        # 需要保留的文件  
        features_pkl = os.path.join(output_dir, f"{name}_features.pkl")  
        masif_pkl = os.path.join(output_dir, f"{name}_features_with_masif.pkl")  
          
        # 检查必要文件是否存在  
        if not os.path.exists(masif_pkl):  
            return False, f"MaSIF PKL file not found: {masif_pkl}"  
              
        # 删除不需要的子目录和文件  
        dirs_to_remove = ['descriptors', 'precomputed', 'surfaces', 'pdbs', 'tmp']  
        for dir_name in dirs_to_remove:  
            dir_path = os.path.join(output_dir, dir_name)  
            if os.path.exists(dir_path):  
                shutil.rmtree(dir_path)  
                print(f"Removed directory: {dir_path}")  
          
        # 删除其他临时文件，但保留必要的PKL文件  
        for file_path in glob.glob(os.path.join(output_dir, "*")):  
            if os.path.isfile(file_path):  
                if file_path not in [features_pkl, masif_pkl]:  
                    os.remove(file_path)  
                    print(f"Removed file: {file_path}")  
          
        print(f"Cleanup completed for {name}")  
        return True, "Success"  
          
    except Exception as e:  
        error_msg = f"Error cleaning up {target_dir}: {e}"  
        print(error_msg)  
        return False, error_msg  
  
def load_pkl_file(pkl_path):  
    """  
    加载单个PKL文件  
    """  
    try:  
        with open(pkl_path, 'rb') as f:  
            data = pickle.load(f)  
        return data  
    except Exception as e:  
        print(f"Error loading {pkl_path}: {e}")  
        return None  
  
def merge_pkl_files(pkl_files, output_file):  
    """  
    合并多个PKL文件  
    """  
    try:  
        all_complexes = []  
          
        print(f"Merging {len(pkl_files)} PKL files...")  
          
        # 使用多进程加载PKL文件  
        with ProcessPoolExecutor(max_workers=24) as executor:  
            results = list(executor.map(load_pkl_file, pkl_files))  
          
        for i, data in enumerate(results):  
            if data is None:  
                print(f"Skipping failed file: {pkl_files[i]}")  
                continue  
                  
            if isinstance(data, list):  
                all_complexes.extend(data)  
            elif isinstance(data, dict):  
                all_complexes.append(data)  
            else:  
                print(f"Unknown data format in {pkl_files[i]}")  
          
        print(f"Total complexes to merge: {len(all_complexes)}")  
          
        # 保存合并后的文件  
        with open(output_file, 'wb') as f:  
            pickle.dump(all_complexes, f)  
              
        print(f"Merged PKL saved to: {output_file}")  
        return True  
          
    except Exception as e:  
        print(f"Error merging PKL files: {e}")  
        traceback.print_exc()  
        return False  
  
def process_csv_and_merge(csv_file, base_output_dir):  
    """  
    处理CSV文件，清理目录并合并PKL文件  
    """  
    try:  
        # 读取CSV文件获取所有复合物信息  
        complexes_info = []  
        with open(csv_file, 'r') as f:  
            lines = f.readlines()  
              
        # 跳过标题行  
        for line in lines[1:]:  
            line = line.strip()  
            if not line:  
                continue  
            parts = line.split(',')  
            if len(parts) >= 3:  
                receptor = parts[0].strip()  
                name = parts[2].strip()  
                target_dir = os.path.dirname(receptor)  
                complexes_info.append((target_dir, name))  
          
        print(f"Found {len(complexes_info)} complexes to process")  
          
        # 记录成功和失败的目录  
        success_dirs = []  
        failed_dirs = []  
          
        # 清理目录并收集PKL文件路径  
        features_pkl_files = []  
        masif_pkl_files = []  
          
        for target_dir, name in complexes_info:  
            print(f"Processing {name}...")  
              
            # 清理目录  
            success, message = cleanup_directory(target_dir, name)  
            if success:  
                success_dirs.append((target_dir, name, message))  
                # 收集PKL文件路径  
                output_dir = os.path.join(target_dir, "output")  
                features_pkl = os.path.join(output_dir, f"{name}_features.pkl")  
                masif_pkl = os.path.join(output_dir, f"{name}_features_with_masif.pkl")  
                  
                if os.path.exists(features_pkl):  
                    features_pkl_files.append(features_pkl)  
                if os.path.exists(masif_pkl):  
                    masif_pkl_files.append(masif_pkl)  
            else:  
                failed_dirs.append((target_dir, name, message))  
          
        # 写入成功和失败的记录文件  
        success_log = os.path.join(base_output_dir, "successful_directories.txt")  
        failed_log = os.path.join(base_output_dir, "failed_directories.txt")  
          
        with open(success_log, 'w') as f:  
            f.write(f"# Successful directories ({len(success_dirs)} total)\n")  
            f.write("# Format: directory_path | complex_name | status\n\n")  
            for target_dir, name, message in success_dirs:  
                f.write(f"{target_dir} | {name} | {message}\n")  
          
        with open(failed_log, 'w') as f:  
            f.write(f"# Failed directories ({len(failed_dirs)} total)\n")  
            f.write("# Format: directory_path | complex_name | error_message\n\n")  
            for target_dir, name, message in failed_dirs:  
                f.write(f"{target_dir} | {name} | {message}\n")  
          
        print(f"Success log saved to: {success_log}")  
        print(f"Failed log saved to: {failed_log}")  
        print(f"Collected {len(features_pkl_files)} features.pkl files")  
        print(f"Collected {len(masif_pkl_files)} features_with_masif.pkl files")  
          
        # 合并PKL文件  
        if features_pkl_files:  
            features_output = os.path.join(base_output_dir, "merged_features.pkl")  
            merge_pkl_files(features_pkl_files, features_output)  
          
        if masif_pkl_files:  
            masif_output = os.path.join(base_output_dir, "merged_features_with_masif.pkl")  
            merge_pkl_files(masif_pkl_files, masif_output)  
          
        # 写入处理总结  
        summary_file = os.path.join(base_output_dir, "processing_summary.txt")  
        with open(summary_file, 'w') as f:  
            f.write(f"Processing Summary\n")  
            f.write(f"==================\n\n")  
            f.write(f"Total complexes processed: {len(complexes_info)}\n")  
            f.write(f"Successful: {len(success_dirs)}\n")  
            f.write(f"Failed: {len(failed_dirs)}\n\n")  
            f.write(f"PKL files merged:\n")  
            f.write(f"- features.pkl: {len(features_pkl_files)} files\n")  
            f.write(f"- features_with_masif.pkl: {len(masif_pkl_files)} files\n\n")  
            f.write(f"Output files:\n")  
            if features_pkl_files:  
                f.write(f"- merged_features.pkl\n")  
            if masif_pkl_files:  
                f.write(f"- merged_features_with_masif.pkl\n")  
            f.write(f"- successful_directories.txt\n")  
            f.write(f"- failed_directories.txt\n")  
          
        print(f"Processing summary saved to: {summary_file}")  
        return True  
          
    except Exception as e:  
        print(f"Error processing CSV and merging: {e}")  
        traceback.print_exc()  
        return False  
  
def main():  
    parser = argparse.ArgumentParser(description='清理临时文件并合并所有PKL文件')  
    parser.add_argument('--csv_file', required=True, help='输入的CSV文件路径')  
    parser.add_argument('--output_dir', required=True, help='合并后PKL文件的输出目录')  
      
    args = parser.parse_args()  
      
    if not os.path.exists(args.csv_file):  
        print(f"Error: CSV file not found: {args.csv_file}")  
        sys.exit(1)  
      
    # 创建输出目录  
    os.makedirs(args.output_dir, exist_ok=True)  
      
    # 处理CSV并合并PKL文件  
    success = process_csv_and_merge(args.csv_file, args.output_dir)  
      
    if not success:  
        print("Processing failed!")  
        sys.exit(1)  
      
    print("All processing completed successfully!")  
  
if __name__ == "__main__":  
    main()