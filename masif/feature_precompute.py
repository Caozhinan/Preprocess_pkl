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
    Step 3: 特征预计算  
      
    参数:  
    - ply_file_path: Step 1和2生成的PLY文件路径  
    - output_dir: 输出目录  
    - masif_app: 应用类型 ('masif_site' 或 'masif_ppi_search')  
      
    返回:  
    - 预计算的特征字典  
    """  
      
    print(f"开始Step 3: 特征预计算...")  
    # print(f"处理PLY文件: {ply_file_path}")  
      
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
      
    # print(f"输出目录: {my_precomp_dir}")  
      
    # 3. 读取表面数据并计算特征  
    # print("读取表面数据并计算patch特征...")  
      
    try:  
        # 调用核心特征计算函数  
        input_feat, rho, theta, mask, neigh_indices, iface_labels, verts = read_data_from_surface(  
            ply_file_path, params  
        )  
          
        # print(f"成功计算特征:")  
        # print(f"  - 顶点数: {len(verts)}")  
        # print(f"  - 输入特征维度: {input_feat.shape}")  
        # # print(f"  - 极坐标 rho: {rho.shape}")  
        # # print(f"  - 极坐标 theta: {theta.shape}")  
        # # print(f"  - 邻居索引: {len(neigh_indices)}")  
          
    except Exception as e:  
        print(f"特征计算失败: {e}")  
        return None  
      
    # 4. 保存预计算的特征  
    print("保存预计算特征...")  
      
    # 保存主要特征数据  
    np.save(os.path.join(my_precomp_dir, 'p1_rho_wrt_center'), rho)  
    np.save(os.path.join(my_precomp_dir, 'p1_theta_wrt_center'), theta)  
    np.save(os.path.join(my_precomp_dir, 'p1_input_feat'), input_feat)  
    np.save(os.path.join(my_precomp_dir, 'p1_mask'), mask)  
    np.save(os.path.join(my_precomp_dir, 'p1_list_indices'), neigh_indices)  
    np.save(os.path.join(my_precomp_dir, 'p1_iface_labels'), iface_labels)  
      
    # 保存坐标信息  
    np.save(os.path.join(my_precomp_dir, 'p1_X.npy'), verts[:, 0])  
    np.save(os.path.join(my_precomp_dir, 'p1_Y.npy'), verts[:, 1])  
    np.save(os.path.join(my_precomp_dir, 'p1_Z.npy'), verts[:, 2])  
      
    # print("特征保存完成!")  
      
    # 5. 构建返回的特征字典  
    features = {  
        'input_features': input_feat,  # 5维特征向量 (shape_index, DDC, hbond, charge, hphob)  
        'rho_coordinates': rho,        # 径向坐标  
        'theta_coordinates': theta,    # 角度坐标  
        'mask': mask,                  # patch掩码  
        'neighbor_indices': neigh_indices,  # 邻居索引  
        'interface_labels': iface_labels,   # 界面标签  
        'vertices': verts,             # 顶点坐标  
        'output_dir': my_precomp_dir,  # 输出目录  
        'num_patches': len(verts),     # patch数量  
        'feature_dim': input_feat.shape[-1]  # 特征维度  
    }  
      
    print(f"Step 3完成! 生成了 {len(verts)} 个patches，每个patch有 {input_feat.shape[-1]} 维特征")  
      
    return features  
  
def load_precomputed_features(precomp_dir, ppi_pair_id='p1'):  
    """  
    加载预计算的特征  
      
    参数:  
    - precomp_dir: 预计算特征目录  
    - ppi_pair_id: 蛋白质对ID (默认 'p1')  
      
    返回:  
    - 特征字典  
    """  
      
    features = {}  
      
    try:  
        features['input_features'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_input_feat.npy'))  
        features['rho_coordinates'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_rho_wrt_center.npy'))  
        features['theta_coordinates'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_theta_wrt_center.npy'))  
        features['mask'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_mask.npy'))  
        features['neighbor_indices'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_list_indices.npy'), allow_pickle=True)  
        features['interface_labels'] = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_iface_labels.npy'))  
          
        # 加载坐标  
        x = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_X.npy'))  
        y = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_Y.npy'))  
        z = np.load(os.path.join(precomp_dir, f'{ppi_pair_id}_Z.npy'))  
        features['vertices'] = np.column_stack([x, y, z])  
          
        # print(f"成功加载预计算特征:")  
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