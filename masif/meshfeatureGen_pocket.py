#!/usr/bin/env python3
import numpy as np
import os
import pymesh
from default_config.masif_opts import masif_opts
from triangulation.computeMSMS import computeMSMS
from triangulation.fixmesh import fix_mesh
from input_output.save_ply import save_ply
from input_output.protonate import protonate
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS
from triangulation.compute_normal import compute_normal
import time
import argparse


def compute_pocket_surface_features(pocket_pdb, output_dir):
    """
    蛋白质口袋的表面特征计算
    直接处理pocket PDB中的所有残基,不进行链提取
    """

    # 1. 设置临时目录
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    pocket_name = os.path.basename(pocket_pdb).replace('.pdb', '')

    # 2. 质子化PDB文件
    protonated_file = os.path.join(tmp_dir, f"{pocket_name}_protonated.pdb")
    protonate(pocket_pdb, protonated_file)

    # 3. 直接计算MSMS表面(不提取链)
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(
        protonated_file,
        protonate=True,   # 已经质子化过了，保持原参数以兼容函数签名
        ligand_code=None  # 不包含配体
    )

    # 4. 计算化学特征
    vertex_hbond = None
    vertex_hphobicity = None

    if masif_opts['use_hbond']:
        vertex_hbond = computeCharges(
            protonated_file.replace('.pdb', ''),
            vertices1,
            names1,
            None,
            None
        )

    if masif_opts['use_hphob']:
        vertex_hphobicity = computeHydrophobicity(names1, None, None)

    # 5. 网格修复和正则化
    mesh = pymesh.form_mesh(vertices1, faces1)
    regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])

    # 6. 计算几何特征
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    # 7. 将特征分配到正则化网格
    if masif_opts['use_hbond'] and vertex_hbond is not None:
        vertex_hbond = assignChargesToNewMesh(
            regular_mesh.vertices, vertices1, vertex_hbond, masif_opts
        )

    if masif_opts['use_hphob'] and vertex_hphobicity is not None:
        vertex_hphobicity = assignChargesToNewMesh(
            regular_mesh.vertices, vertices1, vertex_hphobicity, masif_opts
        )

    # 8. 计算APBS静电特征
    if masif_opts['use_apbs']:
        try:
            vertex_charges = computeAPBS(
                regular_mesh.vertices,
                protonated_file,
                os.path.join(tmp_dir, pocket_name),
                None  # 不需要mol2文件
            )
        except Exception as e:
            print(f"APBS计算失败: {e}")
            vertex_charges = np.zeros(len(regular_mesh.vertices))
    else:
        vertex_charges = np.zeros(len(regular_mesh.vertices))

    # 9. 计算曲率和形状特征
    feature_mesh = pymesh.form_mesh(regular_mesh.vertices, regular_mesh.faces)
    feature_mesh.add_attribute("vertex_mean_curvature")
    feature_mesh.add_attribute("vertex_gaussian_curvature")

    H = feature_mesh.get_attribute("vertex_mean_curvature")
    K = feature_mesh.get_attribute("vertex_gaussian_curvature")

    elem = np.square(H) - K
    elem[elem < 0] = 1e-8
    k1 = H + np.sqrt(elem)
    k2 = H - np.sqrt(elem)

    shape_index = (k1 + k2) / (k1 - k2)
    shape_index = np.arctan(shape_index) * (2 / np.pi)

    # 10. 保存PLY文件
    ply_dir = os.path.join(output_dir, "surfaces")
    os.makedirs(ply_dir, exist_ok=True)

    output_ply = os.path.join(ply_dir, "pocket.ply")
    save_ply(
        output_ply,
        regular_mesh.vertices,
        regular_mesh.faces,
        normals=vertex_normal,
        charges=vertex_charges,
        normalize_charges=True,
        hbond=vertex_hbond,
        hphob=vertex_hphobicity,
        iface=np.zeros(len(regular_mesh.vertices))
    )

    return output_ply


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--pocket_pdb', required=True, help='蛋白质口袋PDB文件')
    parser.add_argument('--output_dir', required=True, help='输出目录')

    args = parser.parse_args()

    result = compute_pocket_surface_features(
        pocket_pdb=args.pocket_pdb,
        output_dir=args.output_dir
    )

    elapsed = time.time() - start_time
    print(f"总耗时: {elapsed:.2f} 秒")