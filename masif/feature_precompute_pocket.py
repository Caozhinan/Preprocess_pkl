import numpy as np
import os
import time
import argparse

from default_config.masif_opts import masif_opts
from masif_modules.read_data_from_surface import read_data_from_surface


def precompute_pocket_features(
    ply_file_path,
    output_dir,
    max_distance=6.0,
    max_shape_size=60
):
    """
    蛋白质口袋的特征预计算

    参数:
    - ply_file_path: 口袋PLY文件路径
    - output_dir: 输出目录
    - max_distance: Patch半径 (默认6.0Å,适合口袋)
    - max_shape_size: 每个patch的最大顶点数 (默认60)
    """

    # 1. 设置自定义参数
    params = {
        'max_distance': max_distance,
        'max_shape_size': max_shape_size,
        'feat_mask': [1.0] * 5
    }

    # 2. 创建输出目录
    pocket_name = os.path.basename(ply_file_path).replace('.ply', '')
    my_precomp_dir = os.path.join(output_dir, 'precomputed', pocket_name + '/')
    os.makedirs(my_precomp_dir, exist_ok=True)

    # 3. 读取表面数据并计算特征
    try:
        input_feat, rho, theta, mask, neigh_indices, iface_labels, verts = \
            read_data_from_surface(ply_file_path, params)
    except Exception as e:
        print(f"特征计算失败: {e}")
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

    return my_precomp_dir


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--ply_file', required=True, help='口袋PLY文件')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    parser.add_argument('--max_distance', type=float, default=6.0,
                        help='Patch半径/Å')
    parser.add_argument('--max_shape_size', type=int, default=60,
                        help='Patch最大顶点数')

    args = parser.parse_args()

    result = precompute_pocket_features(
        ply_file_path=args.ply_file,
        output_dir=args.output_dir,
        max_distance=args.max_distance,
        max_shape_size=args.max_shape_size
    )

    elapsed = time.time() - start_time
    print(f"总耗时: {elapsed:.2f} 秒")