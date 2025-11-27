import numpy as np
import os
import sys
import time
import argparse

# 导入必要的模块
from default_config.masif_opts import masif_opts
from masif_modules.read_data_from_surface import read_data_from_surface
from geometry.compute_polar_coordinates import compute_polar_coordinates


def precompute_ligand_features(
    ply_file_path,
    output_dir,
    max_distance=5.0,
    max_shape_size=60
):
    """
    小分子配体的特征预计算

    参数:
    - ply_file_path: 小分子PLY文件路径 (如 2tpi_ligand.ply)
    - output_dir: 输出目录
    - max_distance: Patch半径 (默认5.0Å,适合小分子)
    - max_shape_size: 每个patch的最大顶点数 (默认60)

    返回:
    - 预计算的特征字典
    """

    # 1. 设置自定义参数
    params = {
        'max_distance': max_distance,
        'max_shape_size': max_shape_size,
        'feat_mask': [1.0] * 5  # 使用所有5维特征
    }

    # 2. 创建输出目录
    ligand_name = os.path.basename(ply_file_path).replace('.ply', '')
    my_precomp_dir = os.path.join(output_dir, 'precomputed', ligand_name + '/')
    os.makedirs(my_precomp_dir, exist_ok=True)

    # 3. 读取表面数据并计算特征
    try:
        input_feat, rho, theta, mask, neigh_indices, iface_labels, verts = \
            read_data_from_surface(ply_file_path, params)
    except Exception as e:
        print(f"特征计算失败: {e}")
        import traceback
        traceback.print_exc()
        return None

    # 4. 保存预计算的特征
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

    np.savez_compressed(
        os.path.join(my_precomp_dir, 'all_features.npz'),
        **features_dict
    )

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

    return features


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(description='小分子配体的特征预计算')
    parser.add_argument('--ply_file', required=True, help='输入PLY文件路径')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    parser.add_argument(
        '--max_distance',
        type=float,
        default=5.0,
        help='Patch半径/Å (默认5.0)'
    )
    parser.add_argument(
        '--max_shape_size',
        type=int,
        default=60,
        help='每个patch的最大顶点数 (默认60)'
    )

    args = parser.parse_args()

    features = precompute_ligand_features(
        ply_file_path=args.ply_file,
        output_dir=args.output_dir,
        max_distance=args.max_distance,
        max_shape_size=args.max_shape_size
    )

    elapsed = time.time() - start_time
    print(f"总耗时: {elapsed:.2f} 秒")

    if features is not None:
        print("小分子特征预计算成功!")
    else:
        print("小分子特征预计算失败!")
        sys.exit(1)