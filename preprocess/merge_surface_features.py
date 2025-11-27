#!/usr/bin/env python3  
"""  
将小分子和蛋白质口袋的表面特征整合到PKL文件中  
过滤掉不必要的字段(pro_name, AA_name)  
"""  
  
import argparse  
import pickle  
import numpy as np  
import os  
import sys  
import traceback  
  
def merge_surface_features_to_pkl(pkl_file_path, ligand_npz_path, pocket_npz_path, output_file=None):  
    """  
    将小分子和蛋白质口袋的p1_input_feat整合到现有PKL文件中  
    同时过滤掉pro_name和AA_name字段  
      
    参数:  
    - pkl_file_path: 原始PKL文件路径 (${name}_features.pkl)  
    - ligand_npz_path: 小分子特征文件路径 (output/precomputed/ligand/all_features.npz)  
    - pocket_npz_path: 蛋白质口袋特征文件路径 (output/precomputed/pocket/all_features.npz)  
    - output_file: 输出PKL文件路径 (默认: ${name}_features_with_masif.pkl)  
      
    返回:  
    - True: 成功  
    - False: 失败  
    """  
    try:  
        # 1. 加载原始PKL文件  
        if not os.path.exists(pkl_file_path):  
            print(f"错误: PKL文件不存在: {pkl_file_path}")  
            return False  
          
        with open(pkl_file_path, 'rb') as f:  
            original_data = pickle.load(f)  
          
        # 2. 加载小分子特征  
        if not os.path.exists(ligand_npz_path):  
            print(f"错误: 小分子NPZ文件不存在: {ligand_npz_path}")  
            return False  
          
        ligand_npz = np.load(ligand_npz_path)  
        if 'p1_input_feat' not in ligand_npz:  
            print(f"错误: 小分子NPZ文件中没有p1_input_feat")  
            return False  
          
        ligand_feat = ligand_npz['p1_input_feat']  
        print(f"小分子特征形状: {ligand_feat.shape}")  
          
        # 3. 加载蛋白质口袋特征  
        if not os.path.exists(pocket_npz_path):  
            print(f"错误: 蛋白质口袋NPZ文件不存在: {pocket_npz_path}")  
            return False  
          
        pocket_npz = np.load(pocket_npz_path)  
        if 'p1_input_feat' not in pocket_npz:  
            print(f"错误: 蛋白质口袋NPZ文件中没有p1_input_feat")  
            return False  
          
        pocket_feat = pocket_npz['p1_input_feat']  
        print(f"蛋白质口袋特征形状: {pocket_feat.shape}")  
          
        # 4. 处理数据结构(列表或字典)  
        if isinstance(original_data, list):  
            if len(original_data) == 0:  
                print("错误: PKL文件为空列表")  
                return False  
              
            # 处理列表中的第一个元素  
            data_dict = original_data[0]  
              
            # 定义需要保留的字段(排除pro_name和AA_name)  
            keep_fields = [  
                'edge_index', 'edge_feat', 'node_feat', 'coords',  
                'smiles', 'rmsd', 'gbscore', 'pk', 'pdbid',  
                'num_node', 'num_edge',  
                'lig_spatial_edge_index', 'lig_spatial_edge_attr',  
                'pro_spatial_edge_index', 'pro_spatial_edge_attr'  
            ]  
              
            # 创建新的字典,只包含需要保留的字段  
            merged_dict = {}  
            for key in keep_fields:  
                if key in data_dict:  
                    merged_dict[key] = data_dict[key]  
              
            # 添加MaSIF特征  
            merged_dict['ligand_masif_feature'] = ligand_feat  
            merged_dict['protein_masif_feature'] = pocket_feat  
              
            # 包装回列表  
            merged_data = [merged_dict]  
              
        elif isinstance(original_data, dict):  
            # 定义需要保留的字段  
            keep_fields = [  
                'edge_index', 'edge_feat', 'node_feat', 'coords',  
                'smiles', 'rmsd', 'gbscore', 'pk', 'pdbid',  
                'num_node', 'num_edge',  
                'lig_spatial_edge_index', 'lig_spatial_edge_attr',  
                'pro_spatial_edge_index', 'pro_spatial_edge_attr'  
            ]  
              
            # 创建新的字典  
            merged_data = {}  
            for key in keep_fields:  
                if key in original_data:  
                    merged_data[key] = original_data[key]  
              
            # 添加MaSIF特征  
            merged_data['ligand_masif_feature'] = ligand_feat  
            merged_data['protein_masif_feature'] = pocket_feat  
              
        else:  
            print(f"错误: 不支持的数据类型: {type(original_data)}")  
            return False  
          
        # 5. 保存合并后的PKL文件  
        output_file = output_file or pkl_file_path.replace('_features.pkl', '_features_with_masif.pkl')  
          
        with open(output_file, 'wb') as f:  
            pickle.dump(merged_data, f, protocol=pickle.HIGHEST_PROTOCOL)  
          
        print(f"合并后的PKL已保存: {output_file}")  
          
        # 6. 验证输出文件  
        if not os.path.exists(output_file):  
            print("错误: 输出文件未成功创建")  
            return False  
          
        # 打印文件大小  
        file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB  
        print(f"输出文件大小: {file_size:.2f} MB")  
          
        return True  
          
    except Exception as e:  
        print(f"合并失败: {e}")  
        traceback.print_exc()  
        return False  
  
def main():  
    parser = argparse.ArgumentParser(description='将表面特征整合到PKL文件(过滤pro_name和AA_name)')  
    parser.add_argument('--pkl_file', required=True, help='原始PKL文件路径')  
    parser.add_argument('--ligand_npz', required=True, help='小分子特征NPZ文件路径')  
    parser.add_argument('--pocket_npz', required=True, help='蛋白质口袋特征NPZ文件路径')  
    parser.add_argument('--output_file', help='输出PKL文件路径(可选)')  
      
    args = parser.parse_args()  
      
    # 验证输入文件  
    if not os.path.exists(args.pkl_file):  
        print(f"错误: PKL文件不存在: {args.pkl_file}")  
        sys.exit(1)  
      
    if not os.path.exists(args.ligand_npz):  
        print(f"错误: 小分子NPZ文件不存在: {args.ligand_npz}")  
        sys.exit(1)  
      
    if not os.path.exists(args.pocket_npz):  
        print(f"错误: 蛋白质口袋NPZ文件不存在: {args.pocket_npz}")  
        sys.exit(1)  
      
    # 执行合并  
    success = merge_surface_features_to_pkl(  
        pkl_file_path=args.pkl_file,  
        ligand_npz_path=args.ligand_npz,  
        pocket_npz_path=args.pocket_npz,  
        output_file=args.output_file  
    )  
      
    if not success:  
        print("合并失败!")  
        sys.exit(1)  
      
    print("合并成功!")  
  
if __name__ == "__main__":  
    main()