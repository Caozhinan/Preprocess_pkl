#!/usr/bin/env python3
import numpy as np
import os
import sys
import argparse
import time
import tempfile
import shutil

# 本地导入
from default_config.masif_opts import masif_opts
from triangulation.computeMSMS import computeMSMS  # 未使用但保留
from triangulation.fixmesh import fix_mesh
from triangulation.ligand_utils import extract_ligand, sdf_to_mol2
import pymesh
from input_output.save_ply import save_ply
from input_output.protonate import protonate  # 未使用但保留
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS  # 未使用但保留
from triangulation.compute_normal import compute_normal
from triangulation.computeAPBS_openbabel import computeAPBS_hybrid, has_unsupported_atoms  # 未使用但保留


def compute_ligand_surface_features(
    sdf_file,
    output_dir,
    density=8.0,
    probe_radius=1.2,
    mesh_res=0.8
):
    """
    专门用于小分子配体的表面特征计算
    """

    tmp_dir = output_dir + "/tmp/" if output_dir else tempfile.gettempdir()
    os.makedirs(tmp_dir, exist_ok=True)

    ligand_name = os.path.basename(sdf_file).replace('.sdf', '')

    # 1. 从 SDF 提取配体
    rdmol = extract_ligand(ligand_sdf_file=sdf_file)
    if rdmol is None:
        print("错误: 无法从SDF文件读取配体")
        return None

    # 2. SDF -> PDB（保留所有氢）
    ligand_pdb = os.path.join(tmp_dir, f"{ligand_name}.pdb")

    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
    if mol is None:
        print("错误: 无法读取SDF分子")
        return None

    writer = Chem.PDBWriter(ligand_pdb)
    writer.write(mol)
    writer.close()

    # 3. 准备 PDB 文件（不额外质子化）
    protonated_file = ligand_pdb

    # 生成 MOL2 文件用于 APBS（即使当前没用到也按原逻辑执行）
    mol2_file = os.path.join(tmp_dir, f"{ligand_name}.mol2")
    if not sdf_to_mol2(sdf_file, mol2_file):
        mol2_file = None

    # 4. 计算 MSMS 表面
    result = computeMSMS_custom(
        protonated_file,
        density=density,
        probe_radius=probe_radius,
        ligand_code=None  # 小分子不需要指定 ligand_code
    )

    if result[0] is None:
        print("错误: MSMS表面生成失败")
        return None

    vertices1, faces1, normals1, names1, areas1 = result

    # 5. 计算化学特征
    vertex_hbond = None
    if masif_opts['use_hbond']:
        vertex_hbond = computeCharges(
            os.path.join(tmp_dir, ligand_name),
            vertices1,
            names1,
            "LIG",
            rdmol
        )

    vertex_hphobicity = None
    if masif_opts['use_hphob']:
        vertex_hphobicity = computeHydrophobicity(names1, "LIG", rdmol)

    # 6. 网格修复和正则化
    mesh = pymesh.form_mesh(vertices1, faces1)
    regular_mesh = fix_mesh(mesh, mesh_res)

    # 7. 计算几何特征（法向）
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    # 8. 将化学特征映射到正则化后网格
    if vertex_hbond is not None:
        vertex_hbond = assignChargesToNewMesh(
            regular_mesh.vertices, vertices1, vertex_hbond, masif_opts
        )

    if vertex_hphobicity is not None:
        vertex_hphobicity = assignChargesToNewMesh(
            regular_mesh.vertices, vertices1, vertex_hphobicity, masif_opts
        )

    # 9. 电荷特征（Gasteiger）
    if masif_opts['use_apbs']:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from sklearn.neighbors import NearestNeighbors

            mol = Chem.SDMolSupplier(sdf_file, removeHs=False)[0]
            if mol is not None:
                AllChem.ComputeGasteigerCharges(mol)
                atom_charges = []
                for atom in mol.GetAtoms():
                    charge = atom.GetDoubleProp('_GasteigerCharge')
                    if np.isnan(charge):
                        charge = 0.0
                    atom_charges.append(charge)

                atom_coords = mol.GetConformer().GetPositions()
                nbrs = NearestNeighbors(n_neighbors=1).fit(atom_coords)
                distances, indices = nbrs.kneighbors(regular_mesh.vertices)

                vertex_charges = np.array(
                    [atom_charges[idx[0]] for idx in indices]
                )
                vertex_charges = np.clip(vertex_charges, -1, 1)
            else:
                vertex_charges = np.zeros(len(regular_mesh.vertices))
        except Exception as e:
            print(f"  电荷计算失败: {e}")
            import traceback
            traceback.print_exc()
            vertex_charges = np.zeros(len(regular_mesh.vertices))
    else:
        vertex_charges = np.zeros(len(regular_mesh.vertices))

    # 10. 曲率与形状特征
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

    # 11. 保存 PLY
    output_ply = os.path.join(tmp_dir, f"{ligand_name}.ply")

    save_ply(
        output_ply,
        regular_mesh.vertices,
        regular_mesh.faces,
        normals=vertex_normal,
        charges=vertex_charges,
        normalize_charges=True,
        hbond=vertex_hbond if vertex_hbond is not None
              else np.zeros(len(regular_mesh.vertices)),
        hphob=vertex_hphobicity if vertex_hphobicity is not None
              else np.zeros(len(regular_mesh.vertices)),
        iface=np.zeros(len(regular_mesh.vertices))
    )

    # 12. 复制到输出目录
    if output_dir:
        ply_dir = os.path.join(output_dir, "surfaces")
        os.makedirs(ply_dir, exist_ok=True)

        final_ply = os.path.join(ply_dir, "ligand.ply")
        shutil.copy(output_ply, final_ply)
    else:
        final_ply = output_ply

    # 13. 构建结果
    result = {
        'vertices': regular_mesh.vertices,
        'faces': regular_mesh.faces,
        'normals': vertex_normal,
        'charges': vertex_charges,
        'hbond': vertex_hbond if vertex_hbond is not None
                 else np.zeros(len(regular_mesh.vertices)),
        'hydrophobicity': vertex_hphobicity if vertex_hphobicity is not None
                          else np.zeros(len(regular_mesh.vertices)),
        'shape_index': shape_index,
        'mean_curvature': H,
        'gaussian_curvature': K,
        'ply_file': final_ply
    }

    # 14. 清理中间文件
    cleanup_files = [
        os.path.join(tmp_dir, f"{ligand_name}.pdb"),
        os.path.join(tmp_dir, f"{ligand_name}.mol2"),
        os.path.join(tmp_dir, f"{ligand_name}.ply"),
        os.path.join(tmp_dir, f"{ligand_name}.csv"),
        os.path.join(tmp_dir, "io.mc"),
    ]

    import glob
    msms_files = glob.glob(os.path.join(tmp_dir, "msms_*.xyzrn")) + \
                 glob.glob(os.path.join(tmp_dir, "msms_*.vert")) + \
                 glob.glob(os.path.join(tmp_dir, "msms_*.face")) + \
                 glob.glob(os.path.join(tmp_dir, "msms_*.area"))
    cleanup_files.extend(msms_files)

    for file_path in cleanup_files:
        try:
            if os.path.exists(file_path):
                os.remove(file_path)
        except Exception:
            # 清理失败不影响主流程，不再打印
            pass

    return result


def computeMSMS_custom(pdb_file, density=8.0, probe_radius=1.2, ligand_code=None):
    import random
    from input_output.read_msms import read_msms
    from triangulation.xyzrn import output_pdb_as_xyzrn
    from default_config.global_vars import msms_bin
    from subprocess import Popen, PIPE

    # 使用 output/*/tmp 目录
    output_dir = os.path.dirname(os.path.dirname(pdb_file))
    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    randnum = random.randint(1, 10000000)
    file_base = os.path.join(tmp_dir, f"msms_{randnum}")
    out_xyzrn = file_base + ".xyzrn"

    ligand_list = ['UNL', 'LIG']
    output_pdb_as_xyzrn(pdb_file, out_xyzrn, keep_hetatms=ligand_list)

    if not os.path.exists(out_xyzrn):
        print(f"错误: XYZRN文件未生成: {out_xyzrn}")
        return None, None, None, None, {}

    with open(out_xyzrn, 'r') as f:
        lines = f.readlines()
        if len(lines) < 2:
            print("错误: XYZRN文件内容不足")
            return None, None, None, None, {}

    args = [
        msms_bin,
        "-density", str(density),
        "-hdensity", str(density),
        "-probe", str(probe_radius),
        "-if", out_xyzrn,
        "-of", file_base,
        "-af", file_base
    ]

    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    if not os.path.exists(file_base + ".vert"):
        print("错误: MSMS未生成 .vert 文件")
        return None, None, None, None, {}

    vertices, faces, normals, names = read_msms(file_base)
    areas = {}

    try:
        ses_file = open(file_base + ".area")
        next(ses_file)
        for line in ses_file:
            fields = line.split()
            areas[fields[3]] = fields[1]
        ses_file.close()
    except Exception:
        pass

    # 临时文件清理在外层统一做
    return vertices, faces, normals, names, areas


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description='生成小分子配体的表面网格和特征'
    )
    parser.add_argument('--sdf_file', required=True, help='输入SDF文件路径')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    parser.add_argument('--density', type=float, default=8.0,
                        help='MSMS表面密度 (默认8.0)')
    parser.add_argument('--probe', type=float, default=1.2,
                        help='探针半径/Å (默认1.2)')
    parser.add_argument('--mesh_res', type=float, default=0.8,
                        help='网格分辨率 (默认0.8)')

    args = parser.parse_args()

    result = compute_ligand_surface_features(
        sdf_file=args.sdf_file,
        output_dir=args.output_dir,
        density=args.density,
        probe_radius=args.probe,
        mesh_res=args.mesh_res
    )

    elapsed = time.time() - start_time
    print(f"总耗时: {elapsed:.2f} 秒")

    if result is not None:
        print("小分子表面特征计算成功!")
    else:
        print("小分子表面特征计算失败!")
        sys.exit(1)