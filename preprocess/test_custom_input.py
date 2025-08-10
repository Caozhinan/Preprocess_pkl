import argparse  # 解析命令行参数
from pathlib import Path  # 跨平台处理文件和路径
import pickle  # 用于保存/加载Python对象
from preprocess import  gen_graph, to_pyg_graph, get_info, GB_score, analyze_plip_interactions  # 预处理和特征生成相关函数
from joblib import Parallel, delayed  # 并行计算工具
from utils import read_mol, obabel_pdb2mol, pymol_pocket  # 读取分子和格式转换工具
import numpy as np  # 数值计算库
from rdkit import Chem, RDLogger  # 化学信息学库及日志管理
import pandas as pd  # 数据处理库
import time
from mol2graph import mol2graph_ligand, mol2graph_protein_from_pdb   # 分子图转换工具
import functools

# 计时装饰器
def timing(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        print(f"[TIMER] {func.__name__} took {elapsed:.3f} seconds")
        return result
    return wrapper

# 给主要步骤打上计时装饰器（假如允许）
read_mol = timing(read_mol)
mol2graph_protein_from_pdb = timing(mol2graph_protein_from_pdb)
get_info = timing(get_info)
mol2graph_ligand = timing(mol2graph_ligand)
GB_score = timing(GB_score)
analyze_plip_interactions = timing(analyze_plip_interactions)
gen_graph = timing(gen_graph)

def parallel_helper(proteinpdb, ligandsdf, name, pk, rmsd, protein_cutoff, pocket_cutoff, spatial_cutoff):  
    RDLogger.DisableLog('rdApp.*')
    t_all = time.perf_counter()

    if not (proteinpdb.is_file() and ligandsdf.is_file()):
        print(f"{proteinpdb} or {ligandsdf} does not exist.")
        return None

    # 生成口袋PDB文件（如果不存在）
    pocketpdb = proteinpdb.parent / (proteinpdb.name.rsplit('.', 1)[0] + '_pocket.pdb')
    if not pocketpdb.is_file():
        t_pocket = time.perf_counter()
        pymol_pocket(proteinpdb, ligandsdf, pocketpdb)
        print(f"[{name}] pymol_pocket: {time.perf_counter()-t_pocket:.3f}s")

    try:
        t0 = time.perf_counter()
        ligand = read_mol(ligandsdf)
        t1 = time.perf_counter()
        pocket_dict = mol2graph_protein_from_pdb(pocketpdb)
        t2 = time.perf_counter()
        proinfo, liginfo = get_info(proteinpdb, ligandsdf)
        t3 = time.perf_counter()
        ligand_dict = mol2graph_ligand(ligand)
        t4 = time.perf_counter()
        res = {
            'lc': ligand_dict['coords'], 'lf': ligand_dict['node_feat'],
            'lei': ligand_dict['edge_index'], 'lea': ligand_dict['edge_feat'],
            'pc': pocket_dict['coords'], 'pf': pocket_dict['node_feat'],
            'pei': pocket_dict['edge_index'], 'pea': pocket_dict['edge_feat'],
            'pdbid': name,
            'ligand_smiles': ligand_dict['smiles'],
            'protein_atom_names': pocket_dict['pro_name'],
            'protein_aa_names': pocket_dict['AA_name']
        }
        t5 = time.perf_counter()
        res['gbscore'] = GB_score(liginfo, proinfo)
        t6 = time.perf_counter()
        plip_interactions = analyze_plip_interactions(str(proteinpdb), str(ligandsdf))
        t7 = time.perf_counter()
        res['plip_interactions'] = plip_interactions
        # 细分步骤打印
        print(f"[{name}] read_mol: {t1-t0:.3f}s, mol2graph_protein_from_pdb: {t2-t1:.3f}s, get_info: {t3-t2:.3f}s, mol2graph_ligand: {t4-t3:.3f}s, GB_score: {t6-t5:.3f}s, plip: {t7-t6:.3f}s")
    except RuntimeError as e:
        print(proteinpdb, pocketpdb, ligandsdf, "Fail on reading molecule")
        return None

    ligand = (res['lc'], res['lf'], res['lei'], res['lea'])
    pocket = (res['pc'], res['pf'], res['pei'], res['pea'])

    try:
        t8 = time.perf_counter()
        raw = gen_graph(
            ligand, pocket, name,
            protein_cutoff=protein_cutoff,
            pocket_cutoff=pocket_cutoff,
            spatial_cutoff=spatial_cutoff,
            protein_file=str(proteinpdb),
            ligand_file=str(ligandsdf),
            plip_interactions=res['plip_interactions']
        )
        t9 = time.perf_counter()
        print(f"[{name}] gen_graph: {t9-t8:.3f}s")
        comp_coord, comp_feat, comp_ei, comp_ea, comp_num_node, comp_num_edge, lig_sei, lig_sea, pro_sei, pro_sea = raw
    except ValueError as e:
        print(f"{name}: Error gen_graph from raw feature {str(e)}")
        return None

    result_dict = {
        'edge_index': comp_ei, 'edge_feat': comp_ea, 'node_feat': comp_feat, 'coords': comp_coord,
        'pro_name': res['protein_atom_names'], 'AA_name': res['protein_aa_names'],
        'smiles': res['ligand_smiles'], 'rmsd': rmsd,
        'gbscore': res['gbscore'],
        'pk': pk, 'pdbid': name, 'num_node': comp_num_node, 'num_edge': comp_num_edge,
        'lig_spatial_edge_index': lig_sei,
        'lig_spatial_edge_attr': lig_sea,
        'pro_spatial_edge_index': pro_sei + len(comp_feat) // 2,  # 蛋白空间边下标偏移
        'pro_spatial_edge_attr': pro_sea
    }
    print(f"[{name}] Total parallel_helper took {time.perf_counter() - t_all:.3f}s")
    return result_dict

if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=Path)  # 输入csv，包含列表
    parser.add_argument('output', type=Path)    # 输出pkl文件路径
    parser.add_argument('--njobs', type=int, default=1)
    parser.add_argument('--protein_cutoff', type=float, default=6.)
    parser.add_argument('--pocket_cutoff', type=float, default=5.)
    parser.add_argument('--spatial_cutoff', type=float, default=5.)
    args = parser.parse_args()

    # 读取csv文件
    filedf = pd.read_csv(args.file_csv)
    receptors = filedf['receptor']
    ligands = filedf['ligand']
    names = filedf['name']
    pks = filedf['pk']
    rmsds = filedf['rmsd']

    graph_dicts = Parallel(n_jobs=args.njobs)(
        delayed(parallel_helper)(
            Path(rec), Path(lig), name, pk, rmsd,
            args.protein_cutoff, args.pocket_cutoff, args.spatial_cutoff
        )
        for rec, lig, name, pk, rmsd in zip(receptors, ligands, names, pks, rmsds)
    )

    # 过滤掉None结果
    graph_dicts = list(filter(None, graph_dicts))

    # 保存为包含所有字典的列表
    pickle.dump(graph_dicts, open(args.output, 'wb'))

    elapsed = time.time() - start_time
    print(f"[custom_input.py] Finished in {elapsed:.2f} seconds")