#!/usr/bin/env python3  
"""  
整合MaSIF指纹到现有PKL文件（仅标准方向指纹）  
"""  
  
import argparse  
import pickle  
import numpy as np  
import os  
import sys  
import traceback  
  
def merge_masif_to_existing_pkl(pkl_file_path, descriptors_dir, output_file=None):  
    """  
    将MaSIF指纹整合到现有的PKL文件中（仅保留标准方向指纹）  
    """  
    try:  
        # 1. 加载现有的PKL文件  
        if not os.path.exists(pkl_file_path):  
            return False  
        with open(pkl_file_path, 'rb') as f:  
            existing_data = pickle.load(f)  
  
        # 2. 检查描述符目录存在性  
        if not os.path.exists(descriptors_dir):  
            return False  
  
        # 3. 初始化MaSIF特征字典  
        masif_features = {}  
  
        # 4. 加载描述符文件（仅标准方向指纹）  
        file_path = os.path.join(descriptors_dir, 'p1_desc_straight.npy')  
        if os.path.exists(file_path):  
            masif_features['masif_desc_straight'] = np.load(file_path)  
        else:  
            return False  
  
        # 5. 整合特征到现有数据中  
        if isinstance(existing_data, dict):  
            existing_data.update(masif_features)  
        elif isinstance(existing_data, list) and len(existing_data) > 0:  
            if isinstance(existing_data[0], dict):  
                existing_data[0].update(masif_features)  
            else:  
                return False  
        else:  
            return False  
  
        # 6. 保存整合后的数据  
        if output_file is None:  
            base_name = os.path.splitext(pkl_file_path)[0]  
            output_file = f"{base_name}_with_masif.pkl"  
  
        with open(output_file, 'wb') as f:  
            pickle.dump(existing_data, f)  
  
        if not os.path.exists(output_file):  
            return False  
  
        return True  
  
    except Exception:  
        traceback.print_exc()  
        return False  
  
def main():  
    parser = argparse.ArgumentParser(description='将MaSIF指纹整合到现有PKL文件（仅标准方向指纹）')  
    parser.add_argument('--pkl_file', required=True, help='现有的PKL文件路径')  
    parser.add_argument('--descriptors_dir', required=True, help='描述符目录')  
    parser.add_argument('--output_file', help='输出PKL文件路径（可选）')  
  
    args = parser.parse_args()  
  
    if not os.path.exists(args.pkl_file):  
        sys.exit(1)  
    if not os.path.exists(args.descriptors_dir):  
        sys.exit(1)  
  
    success = merge_masif_to_existing_pkl(  
        pkl_file_path=args.pkl_file,  
        descriptors_dir=args.descriptors_dir,  
        output_file=args.output_file  
    )  
  
    if not success:  
        sys.exit(1)  
  
if __name__ == "__main__":  
    main()