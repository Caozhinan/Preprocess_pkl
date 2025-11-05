#!/usr/bin/python  
import numpy as np  
import os  
import Bio  
import shutil  
from Bio.PDB import *   
import sys  
import importlib  
from IPython.core.debugger import set_trace  
import argparse
# 本地导入  
from default_config.masif_opts import masif_opts  
from triangulation.computeMSMS import computeMSMS  
from triangulation.fixmesh import fix_mesh  
from triangulation.ligand_utils import extract_ligand, sdf_to_mol2
import pymesh  
from input_output.extractPDB import extractPDB  
from input_output.save_ply import save_ply  
from input_output.read_ply import read_ply  
from input_output.protonate import protonate , fix_ligand_atom_names
from triangulation.computeHydrophobicity import computeHydrophobicity  
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh  
from triangulation.computeAPBS import computeAPBS  
from triangulation.compute_normal import compute_normal  
from sklearn.neighbors import KDTree  
import time

def compute_protein_ligand_surface_features(pdb_file, chain_id, ligand_code=None, ligand_chain=None, sdf_file=None, mol2_patch=None, output_dir=None):  
    """  
    完整的蛋白质-配体表面特征计算流程  
      
    参数:  
    - pdb_file: PDB文件路径  
    - chain_id: 蛋白质链ID  
    - ligand_code: 配体三字母代码 (可选)  
    - ligand_chain: 配体所在链 (可选)  
    - sdf_file: SDF模板文件路径 (可选)  
    - mol2_patch: MOL2补丁文件路径 (可选)  
    - output_dir: 输出目录  
      
    返回:  
    - 增强的特征字典，包含几何、化学和形状特征  
    """  
      
    print("开始蛋白质-配体表面特征计算...")  
      
    # ========== Step 1: PDB提取和表面三角化 ==========  
      
    # 1. 设置临时目录  
    tmp_dir = output_dir + "/tmp/" if output_dir else masif_opts['tmp_dir']  
    os.makedirs(tmp_dir, exist_ok=True)  
      
    # 2. 质子化PDB文件  
    pdb_id = os.path.basename(pdb_file).replace('.pdb', '')  
    protonated_file = tmp_dir + "/" + pdb_id + "_protonated.pdb"  
    # print("质子化PDB文件...")  
    protonate(pdb_file, protonated_file)  
    fix_ligand_atom_names(protonated_file) 
    # 3. 提取指定链  
    out_filename = tmp_dir + "/" + pdb_id + "_" + chain_id  
    print(f"提取链 {chain_id}...")  
    extractPDB(protonated_file, out_filename + ".pdb", chain_id, ligand_code, ligand_chain)  
      
    # 4. 计算MSMS表面  
    # print("计算MSMS分子表面...")  
    vertices1, faces1, normals1, names1, areas1 = computeMSMS(  
        out_filename + ".pdb",   
        protonate=True,   
        ligand_code=ligand_code  
    )  
      
    # 5. 处理配体（如果存在）  
    rdmol = None  
    mol2_file = None  
    if ligand_code is not None and sdf_file is not None:  
        print(f"处理配体 {ligand_code}...")  

        # 检查 SDF 文件是否存在  
        if not os.path.exists(sdf_file):  
            print(f"警告: SDF 文件不存在: {sdf_file}")  
            # print("跳过配体处理...")  
            rdmol = None  
        else:  
            # 调用修改后的 extract_ligand 函数  
            rdmol = extract_ligand(ligand_sdf_file=sdf_file)  

            if rdmol is not None:  
                # 生成 MOL2 文件用于 APBS  
                mol2_file = os.path.join(tmp_dir, "{}_{}.mol2".format(ligand_code, ligand_chain))  

                # 从SDF转换为MOL2  
                if sdf_to_mol2(sdf_file, mol2_file):  
                    print(f"成功生成 MOL2 文件: {mol2_file}")  
                else:  
                    print(f"MOL2 文件生成失败")  
                    mol2_file = None
      
    # ========== Step 2: 表面特征计算 ==========  
      
    # 6. 计算化学特征  
    # print("计算化学特征...")  
      
    vertex_hbond = None  
    vertex_hphobicity = None  
    vertex_charges = None  
      
    # 氢键特征  
    if masif_opts['use_hbond']:  
        # print("  - 计算氢键特征...")  
        vertex_hbond = computeCharges(out_filename, vertices1, names1, ligand_code, rdmol)  
      
    # 疏水性特征  
    if masif_opts['use_hphob']:  
        # print("  - 计算疏水性特征...")  
        vertex_hphobicity = computeHydrophobicity(names1, ligand_code, rdmol)  
      
    # 7. 网格修复和正则化  
    # print("修复和正则化网格...")  
    mesh = pymesh.form_mesh(vertices1, faces1)  
    regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])  
    # print("网格修复完成!")  
      
    # 8. 计算几何特征  
    # print("计算几何特征...")  
      
    # 计算表面法向量  
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)  
      
    # 9. 将特征分配到正则化网格  
    if masif_opts['use_hbond'] and vertex_hbond is not None:  
        vertex_hbond = assignChargesToNewMesh(  
            regular_mesh.vertices, vertices1, vertex_hbond, masif_opts  
        )  
      
    if masif_opts['use_hphob'] and vertex_hphobicity is not None:  
        vertex_hphobicity = assignChargesToNewMesh(  
            regular_mesh.vertices, vertices1, vertex_hphobicity, masif_opts  
        )  
      
    # 10. 计算APBS静电特征（如果启用）  
    if masif_opts['use_apbs']:  
        # print("计算APBS静电特征...")  
        vertex_charges = computeAPBS(  
            regular_mesh.vertices,   
            out_filename + ".pdb",   
            out_filename,   
            mol2_file  
        )  
        print("APBS计算完成!")  
      
    # 11. 计算界面特征（可选）  
    iface = np.zeros(len(regular_mesh.vertices))  
    if 'compute_iface' in masif_opts and masif_opts['compute_iface']:  
        # print("计算界面特征...")  
        # 计算整个复合物的表面  
        v3, f3, _, _, _ = computeMSMS(protonated_file, protonate=True, ligand_code=ligand_code)  
        full_mesh = pymesh.form_mesh(v3, f3)  
          
        # 找到界面顶点  
        kdt = KDTree(full_mesh.vertices)  
        d, r = kdt.query(regular_mesh.vertices)  
        d = np.square(d)  
        iface_v = np.where(d >= 2.0)[0]  
        iface[iface_v] = 1.0  
      
    # 12. 计算曲率和形状特征  
    # print("计算曲率和形状特征...")  
      
    # 创建PyMesh对象以计算曲率  
    feature_mesh = pymesh.form_mesh(regular_mesh.vertices, regular_mesh.faces)  
    feature_mesh.add_attribute("vertex_mean_curvature")  
    feature_mesh.add_attribute("vertex_gaussian_curvature")  
      
    H = feature_mesh.get_attribute("vertex_mean_curvature")  
    K = feature_mesh.get_attribute("vertex_gaussian_curvature")  
      
    # 计算主曲率  
    elem = np.square(H) - K  
    elem[elem < 0] = 1e-8  # 避免负值  
    k1 = H + np.sqrt(elem)  
    k2 = H - np.sqrt(elem)  
      
    # 计算形状指数  
    shape_index = (k1 + k2) / (k1 - k2)  
    shape_index = np.arctan(shape_index) * (2 / np.pi)  
      
    # 13. 保存PLY文件  
    output_ply = out_filename + ".ply"  
    save_ply(  
        output_ply,  
        regular_mesh.vertices,  
        regular_mesh.faces,  
        normals=vertex_normal,  
        charges=vertex_charges,  
        normalize_charges=True,  
        hbond=vertex_hbond,  
        hphob=vertex_hphobicity,  
        iface=iface  
    )  
      
    # 14. 复制到输出目录  
    if output_dir:  
        ply_dir = os.path.join(output_dir, "surfaces")  
        pdb_dir = os.path.join(output_dir, "pdbs")  
        os.makedirs(ply_dir, exist_ok=True)  
        os.makedirs(pdb_dir, exist_ok=True)  
          
        shutil.copy(output_ply, ply_dir)  
        shutil.copy(out_filename + '.pdb', pdb_dir)  
      
    # 15. 构建增强特征字典  
    # print("构建增强特征字典...")  
      
    enhanced_features = {  
        # 基础几何特征  
        'vertices': regular_mesh.vertices,  
        'faces': regular_mesh.faces,  
        'normals': vertex_normal,  
          
        # 化学特征  
        'charges': vertex_charges if vertex_charges is not None else np.zeros(len(regular_mesh.vertices)),  
        'hbond': vertex_hbond if vertex_hbond is not None else np.zeros(len(regular_mesh.vertices)),  
        'hydrophobicity': vertex_hphobicity if vertex_hphobicity is not None else np.zeros(len(regular_mesh.vertices)),  
          
        # 形状特征  
        'shape_index': shape_index,  
        'mean_curvature': H,  
        'gaussian_curvature': K,  
        'principal_curvature_k1': k1,  
        'principal_curvature_k2': k2,  
          
        # 界面特征  
        'interface': iface,  
          
        # 文件路径  
        'ply_file': output_ply,  
        'pdb_file': out_filename + '.pdb',  
          
        # 配体信息  
        'ligand_mol': rdmol,  
        'ligand_code': ligand_code,  
        'ligand_chain': ligand_chain  
    }  
      
    # print("表面特征计算完成!")  
    print(f"顶点数: {len(enhanced_features['vertices'])}")  
    print(f"面数: {len(enhanced_features['faces'])}")  
    print(f"特征维度: {len([k for k in enhanced_features.keys() if isinstance(enhanced_features[k], np.ndarray) and len(enhanced_features[k]) == len(enhanced_features['vertices'])])}")  
      
    return enhanced_features  
  
# 使用示例  
if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_file', required=True)
    parser.add_argument('--chain_id', required=True)
    parser.add_argument('--ligand_code', default=None)
    parser.add_argument('--ligand_chain', default=None)
    parser.add_argument('--sdf_file', default=None)
    parser.add_argument('--output_dir', required=True)

    args = parser.parse_args()

    result = compute_protein_ligand_surface_features(
        pdb_file=args.pdb_file,
        chain_id=args.chain_id,
        ligand_code=args.ligand_code,
        ligand_chain=args.ligand_chain,
        sdf_file=args.sdf_file,
        output_dir=args.output_dir
    )

    elapsed = time.time() - start_time
    print(f"[meshfeatureGen.py] Finished in {elapsed:.2f} seconds")