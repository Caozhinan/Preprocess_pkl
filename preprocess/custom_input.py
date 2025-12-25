#!/usr/bin/env python3
import argparse
from pathlib import Path
import pickle
from preprocess import gen_graph, get_info, GB_score, analyze_plip_interactions
from joblib import Parallel, delayed
from utils import read_mol, pymol_pocket
import numpy as np
from rdkit import RDLogger
from rdkit import Chem
import pandas as pd
from mol2graph import mol2graph_ligand, mol2graph_protein_from_pdb
import time

def parallel_helper(proteinpdb, ligand_file, name, pk, rmsd, protein_cutoff, pocket_cutoff, spatial_cutoff):
    """
    处理单个蛋白质-配体复合物的并行辅助函数
    支持ligand.pdb或ligand.sdf输入
    """
    RDLogger.DisableLog('rdApp.*')
    
    if not (proteinpdb.is_file() and ligand_file.is_file()):
        print(f"{proteinpdb} or {ligand_file} does not exist.")
        return None
    
    # 生成口袋PDB文件（如果不存在）
<<<<<<< HEAD
    pocketpdb = proteinpdb.parent / (proteinpdb.name.rsplit('.', 1)[0] + '_pocket.pdb')
=======
    pocketpdb = proteinpdb.parent / (name + '_protein_pocket.pdb')
>>>>>>> 4f9f0f7 (ppi_to_be_done)
    if not pocketpdb.is_file():
        pymol_pocket(proteinpdb, ligand_file, pocketpdb)
    
    try:
        # 根据文件扩展名选择处理方式
        ligand_suffix = ligand_file.suffix.lower()
        
        if ligand_suffix == '.pdb':
            # 使用mol2graph_protein_from_pdb处理PDB格式配体
            print(f"检测到PDB格式配体: {ligand_file}")
            ligand_dict = mol2graph_protein_from_pdb(ligand_file)
            
            # 尝试从PDB生成SMILES(可选,可能失败)
            try:
                ligand_mol = Chem.MolFromPDBFile(str(ligand_file), removeHs=False)
                if ligand_mol:
                    ligand_dict['smiles'] = Chem.MolToSmiles(ligand_mol)
                else:
                    ligand_dict['smiles'] = 'N/A'
            except:
                ligand_dict['smiles'] = 'N/A'
                
        elif ligand_suffix in ['.sdf', '.mol', '.mol2']:
            # 使用mol2graph_ligand处理SDF格式配体(原有逻辑)
            print(f"检测到SDF格式配体: {ligand_file}")
            ligand = read_mol(ligand_file)
            ligand_dict = mol2graph_ligand(ligand)
        else:
            raise ValueError(f"不支持的配体文件格式: {ligand_suffix}")
        
        # 使用口袋PDB而不是完整蛋白质
        pocket_dict = mol2graph_protein_from_pdb(pocketpdb)
        
        proinfo, liginfo = get_info(proteinpdb, ligand_file)
        
        res = {
            'lc': ligand_dict['coords'], 'lf': ligand_dict['node_feat'],
            'lei': ligand_dict['edge_index'], 'lea': ligand_dict['edge_feat'],
            'pc': pocket_dict['coords'], 'pf': pocket_dict['node_feat'],
            'pei': pocket_dict['edge_index'], 'pea': pocket_dict['edge_feat'],
            'pdbid': name,
            'ligand_smiles': ligand_dict.get('smiles', 'N/A'),
            'protein_atom_names': pocket_dict['pro_name'],
            'protein_aa_names': pocket_dict['AA_name']
        }
        
        res['gbscore'] = GB_score(liginfo, proinfo)
        
        # PLIP相互作用分析
        plip_interactions = analyze_plip_interactions(str(proteinpdb), str(ligand_file))
        if plip_interactions:
            print(f"发现的结合位点数: {len(plip_interactions)}")
        res['plip_interactions'] = plip_interactions
        
    except RuntimeError as e:
        print(f"{proteinpdb}, {pocketpdb}, {ligand_file}: 读取分子失败 - {e}")
        return None
    except ValueError as e:
        print(f"{name}: {e}")
        return None
    
    ligand = (res['lc'], res['lf'], res['lei'], res['lea'])
    pocket = (res['pc'], res['pf'], res['pei'], res['pea'])
    
    try:
        raw = gen_graph(
            ligand, pocket, name,
            protein_cutoff=protein_cutoff,
            pocket_cutoff=pocket_cutoff,
            spatial_cutoff=spatial_cutoff,
            protein_file=str(proteinpdb),
            ligand_file=str(ligand_file),
            plip_interactions=res['plip_interactions'],
            ligand_dict=ligand_dict,
            pocket_dict=pocket_dict
        )
        
        comp_coord, comp_feat, comp_ei, comp_ea, comp_num_node, comp_num_edge, lig_sei, lig_sea, pro_sei, pro_sea = raw
    except ValueError as e:
        print(f"{name}: gen_graph错误 - {str(e)}")
        return None
    
    # 返回字典格式
    result_dict = {
        'edge_index': comp_ei, 'edge_feat': comp_ea, 'node_feat': comp_feat, 'coords': comp_coord,
        'pro_name': res['protein_atom_names'], 'AA_name': res['protein_aa_names'],
        'smiles': res['ligand_smiles'], 'rmsd': rmsd,
        'gbscore': res['gbscore'],
        'pk': pk, 'pdbid': name, 'num_node': comp_num_node, 'num_edge': comp_num_edge,
        'lig_spatial_edge_index': lig_sei,
        'lig_spatial_edge_attr': lig_sea,
        'pro_spatial_edge_index': pro_sei + len(comp_feat) // 2,
        'pro_spatial_edge_attr': pro_sea
    }
    
    return result_dict

if __name__ == "__main__":
    start_time = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=Path, help='输入CSV文件')
    parser.add_argument('output', type=Path, help='输出PKL文件路径')
    parser.add_argument('--njobs', type=int, default=1)
    parser.add_argument('--protein_cutoff', type=float, default=6.)
    parser.add_argument('--pocket_cutoff', type=float, default=5.)
    parser.add_argument('--spatial_cutoff', type=float, default=5.)
    args = parser.parse_args()
    
    # 读取CSV文件
    filedf = pd.read_csv(args.file_csv)
    receptors = filedf['receptor']
    ligands = filedf['ligand']
    names = filedf['name']
    pks = filedf['pk']
    rmsds = filedf['rmsd']
    
    # 并行处理
    graph_dicts = Parallel(n_jobs=args.njobs)(
        delayed(parallel_helper)(
            Path(rec), Path(lig), name, pk, rmsd,
            args.protein_cutoff, args.pocket_cutoff, args.spatial_cutoff
        )
        for rec, lig, name, pk, rmsd in zip(receptors, ligands, names, pks, rmsds)
    )
    
    # 过滤掉None结果
    graph_dicts = list(filter(None, graph_dicts))
    
    # 创建输出目录（如果不存在）
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # 保存为PKL文件
    with open(args.output, 'wb') as f:
        pickle.dump(graph_dicts, f)
    
    elapsed = time.time() - start_time
    print(f"[custom_input.py] 完成,耗时 {elapsed:.2f} 秒")
    print(f"已保存 {len(graph_dicts)} 个复合物到 {args.output}")