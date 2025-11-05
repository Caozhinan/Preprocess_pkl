#!/usr/bin/python  
import numpy as np  
import os  
import sys  
import time  
import pymesh  
from sklearn.neighbors import KDTree  
  
# 导入必要的模块  
from default_config.masif_opts import masif_opts  
from masif_modules.read_data_from_surface import read_data_from_surface, compute_shape_complementarity  
from geometry.compute_polar_coordinates import compute_polar_coordinates  
  
def precompute_surface_features(ply_file_path, output_dir, masif_app='masif_site'):    
    """    
    Step 3: 特征预计算（优化版本 - 使用压缩打包）  
        
    参数:    
    - ply_file_path: Step 1和2生成的PLY文件路径    
    - output_dir: 输出目录    
    - masif_app: 应用类型 ('masif_site' 或 'masif_ppi_search')    
        
    返回:    
    - 预计算的特征字典    
    """    
        
    print(f"开始Step 3: 特征预计算...")    
        
    # 1. 设置参数    
    if masif_app == 'masif_ppi_search':     
        params = masif_opts['ppi_search']    
    elif masif_app == 'masif_site':    
        params = masif_opts['site']    
        params['ply_chain_dir'] = masif_opts['ply_chain_dir']    
    else:    
        raise ValueError(f"不支持的应用类型: {masif_app}")    
        
    # 2. 创建输出目录    
    ppi_pair_id = os.path.basename(ply_file_path).replace('.ply', '')    
    my_precomp_dir = os.path.join(output_dir, 'precomputed', ppi_pair_id + '/')    
    os.makedirs(my_precomp_dir, exist_ok=True)    
        
    # 3. 读取表面数据并计算特征    
    try:    
        # 调用核心特征计算函数    
        input_feat, rho, theta, mask, neigh_indices, iface_labels, verts = read_data_from_surface(    
            ply_file_path, params    
        )    
            
    except Exception as e:    
        print(f"特征计算失败: {e}")    
        return None    
        
    # 4. 保存预计算的特征（优化版本 - 使用压缩打包）  
    # print("保存预计算特征...")    
        
    # 将所有数据打包到一个压缩文件中  
    features_dict = {  
        'p1_rho_wrt_center': rho,  
        'p1_theta_wrt_center': theta,  
        'p1_input_feat': input_feat,  
        'p1_mask': mask,  
        'p1_list_indices': neigh_indices,  
        'p1_iface_labels': iface_labels,  
        'p1_X': verts[:, 0],  
        'p1_Y': verts[:, 1],  
        'p1_Z': verts[:, 2]  
    }  
  
    # 使用压缩格式保存所有特征  
    np.savez_compressed(os.path.join(my_precomp_dir, 'all_features.npz'), **features_dict)  
        
    # 5. 构建返回的特征字典    
    features = {    
        'input_features': input_feat,  
        'rho_coordinates': rho,  
        'theta_coordinates': theta,  
        'mask': mask,  
        'neighbor_indices': neigh_indices,  
        'interface_labels': iface_labels,  
        'vertices': verts,  
        'output_dir': my_precomp_dir,  
        'num_patches': len(verts),  
        'feature_dim': input_feat.shape[-1]  
    }    
        
    print(f"Step 3完成! 生成了 {len(verts)} 个patches，每个patch有 {input_feat.shape[-1]} 维特征")    
        
    return features  
  
def load_precomputed_features(precomp_dir, ppi_pair_id='p1'):    
    """    
    加载预计算的特征（优化版本 - 从压缩文件加载）  
        
    参数:    
    - precomp_dir: 预计算特征目录    
    - ppi_pair_id: 蛋白质对ID (默认 'p1')    
        
    返回:    
    - 特征字典    
    """    
        
    features = {}    
        
    try:    
        # 从压缩文件加载所有特征  
        npz_file = np.load(os.path.join(precomp_dir, 'all_features.npz'))  
          
        features['input_features'] = npz_file['p1_input_feat']  
        features['rho_coordinates'] = npz_file['p1_rho_wrt_center']  
        features['theta_coordinates'] = npz_file['p1_theta_wrt_center']  
        features['mask'] = npz_file['p1_mask']  
        features['neighbor_indices'] = npz_file['p1_list_indices']  
        features['interface_labels'] = npz_file['p1_iface_labels']  
          
        # 重构坐标  
        x = npz_file['p1_X']  
        y = npz_file['p1_Y']  
        z = npz_file['p1_Z']  
        features['vertices'] = np.column_stack([x, y, z])  
          
        npz_file.close()  # 释放文件句柄  
          
        print(f"  - 输入特征: {features['input_features'].shape}")  
        print(f"  - 顶点数: {len(features['vertices'])}")  
          
        return features  
          
    except Exception as e:  
        print(f"加载预计算特征失败: {e}")  
        return None  
  
# 使用示例  
if __name__ == "__main__":
    import argparse

    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--ply_file', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--masif_app', default='masif_site')

    args = parser.parse_args()

    features = precompute_surface_features(
        ply_file_path=args.ply_file,
        output_dir=args.output_dir,
        masif_app=args.masif_app
    )

    elapsed = time.time() - start_time
    print(f"[feature_precompute.py] Finished in {elapsed:.2f} seconds")