#!/usr/bin/env python3    
"""    
整合预处理特征和MaSIF指纹到现有PKL文件    
"""    
  
import argparse    
import pickle    
import numpy as np    
import os    
import sys    
  
def merge_masif_to_existing_pkl(pkl_file, output_dir, name, pk, rmsd, chain_id="A", output_file=None):    
    """    
    将MaSIF指纹整合到现有的PKL文件中    
        
    参数:    
    - pkl_file: 现有的PKL文件路径（custom_input.py生成的）    
    - output_dir: MaSIF输出目录（包含precomputed和descriptors子目录）    
    - name: 复合物名称    
    - pk: pK值    
    - rmsd: RMSD值    
    - chain_id: 蛋白质链标识符，默认为'A'  
    - output_file: 输出PKL文件路径（可选）    
        
    返回:    
    - 成功返回True，失败返回False    
    """    
        
    try:    
        # 加载现有的PKL文件    
        # print(f"加载现有PKL文件: {pkl_file}")    
        with open(pkl_file, 'rb') as f:    
            existing_data = pickle.load(f)    
            
        # print(f"现有PKL文件包含 {len(existing_data)} 个复合物")    
            
        # 设置MaSIF输出路径（动态构建）  
        precomputed_dir = os.path.join(output_dir, 'precomputed', f'complex_{chain_id}')    
        desc_dir = os.path.join(output_dir, 'descriptors', f'complex_{chain_id}')    
            
        # 检查MaSIF输出目录是否存在    
        if not os.path.exists(precomputed_dir):    
            print(f"错误: 预计算特征目录不存在: {precomputed_dir}")    
            return False    
                
        if not os.path.exists(desc_dir):    
            print(f"错误: 指纹目录不存在: {desc_dir}")    
            return False    
            
        # 加载MaSIF特征    
        # print("加载MaSIF特征...")    
        masif_features = {}    
            
        # 预计算特征    
        precomputed_files = [    
            'p1_input_feat.npy',    
            'p1_rho_wrt_center.npy',     
            'p1_theta_wrt_center.npy',    
            'p1_mask.npy'    
        ]    
            
        for file_name in precomputed_files:    
            file_path = os.path.join(precomputed_dir, file_name)    
            if not os.path.exists(file_path):    
                print(f"错误: 缺少预计算特征文件: {file_path}")    
                return False    
                
            key_name = 'masif_' + file_name.replace('p1_', '').replace('.npy', '')    
            masif_features[key_name] = np.load(file_path)    
            # print(f"  加载 {key_name}: {masif_features[key_name].shape}")    
            
        # MaSIF指纹    
        desc_files = [    
            'p1_desc_straight.npy',    
            'p1_desc_flipped.npy'    
        ]    
            
        for file_name in desc_files:    
            file_path = os.path.join(desc_dir, file_name)    
            if not os.path.exists(file_path):    
                print(f"错误: 缺少指纹文件: {file_path}")    
                return False    
                
            key_name = file_name.replace('p1_', 'masif_').replace('.npy', '')    
            masif_features[key_name] = np.load(file_path)    
            # print(f"  加载 {key_name}: {masif_features[key_name].shape}")    
            
        # 找到对应的复合物并添加MaSIF特征    
        target_found = False    
        for i, graph in enumerate(existing_data):    
            # 根据名称匹配（假设PKL中有name或pdbid字段）    
            graph_name = graph.get('name', graph.get('pdbid', ''))    
            if graph_name == name:    
                print(f"找到目标复合物: {name}")    
                # 添加MaSIF特征到现有图数据    
                for key, value in masif_features.items():    
                    graph[key] = value    
                target_found = True    
                break    
            
        if not target_found:    
            print(f"警告: 在PKL文件中未找到名为 {name} 的复合物")    
            # print("可用的复合物名称:")    
            for graph in existing_data:    
                print(f"  - {graph.get('name', graph.get('pdbid', 'Unknown'))}")    
            return False    
            
        # 保存更新后的PKL文件    
        if output_file is None:    
            output_file = pkl_file.replace('.pkl', '_with_masif.pkl')    
            
        print(f"保存整合后的PKL文件到: {output_file}")    
        with open(output_file, 'wb') as f:    
            pickle.dump(existing_data, f)    
            
        # print("整合完成!")    
        # print(f"MaSIF特征总结:")    
        # for key, value in masif_features.items():    
        #     print(f"  - {key}: {value.shape}")    
            
        return True    
            
    except Exception as e:    
        print(f"整合过程中出现错误: {e}")    
        return False    
  
def main():
    import time
    start_time = time.time()

    parser = argparse.ArgumentParser(description='将MaSIF指纹整合到现有PKL文件')
    parser.add_argument('--pkl_file', required=True, help='现有的PKL文件路径')
    parser.add_argument('--output_dir', required=True, help='MaSIF输出目录')
    parser.add_argument('--name', required=True, help='复合物名称')
    parser.add_argument('--pk', required=True, type=float, help='pK值')
    parser.add_argument('--rmsd', required=True, type=float, help='RMSD值')
    parser.add_argument('--chain_id', default='A', help='蛋白质链标识符')
    parser.add_argument('--output_file', help='输出PKL文件路径（可选）')

    args = parser.parse_args()

    # 检查文件是否存在
    if not os.path.exists(args.pkl_file):
        print(f"错误: PKL文件不存在: {args.pkl_file}")
        sys.exit(1)

    if not os.path.exists(args.output_dir):
        print(f"错误: 输出目录不存在: {args.output_dir}")
        sys.exit(1)

    # 执行整合
    success = merge_masif_to_existing_pkl(
        pkl_file=args.pkl_file,
        output_dir=args.output_dir,
        name=args.name,
        pk=args.pk,
        rmsd=args.rmsd,
        chain_id=args.chain_id,
        output_file=args.output_file
    )

    if not success:
        print("整合失败!")
        sys.exit(1)

    print("整合成功!")
    elapsed = time.time() - start_time
    print(f"[merge_pkl.py] Finished in {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()