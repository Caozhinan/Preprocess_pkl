import argparse  # 解析命令行参数
from pathlib import Path  # 跨平台处理文件和路径
# from tqdm import tqdm  # 进度条，可显示任务进度
import pickle  # 用于保存/加载Python对象
from preprocess import  gen_graph, get_info, GB_score, analyze_plip_interactions  # 预处理和特征生成相关函数
from joblib import Parallel, delayed  # 并行计算工具
from utils import read_mol, obabel_pdb2mol, pymol_pocket  # 读取分子和格式转换工具
import numpy as np  # 数值计算库
from rdkit import Chem, RDLogger  # 化学信息学库及日志管理
import tempfile  # 临时文件/目录工具（未使用）
import pandas as pd  # 数据处理库
import os  # 操作系统接口（未使用）
from mol2graph import mol2graph_ligand, mol2graph_protein_from_pdb   # 分子图转换工具
import time

# def process_one(proteinpdb: Path, ligandsdf: Path, name: str, pk: float, rmsd: float, protein_cutoff, pocket_cutoff, spatial_cutoff):  
#     RDLogger.DisableLog('rdApp.*')  
  
#     if not (proteinpdb.is_file() and ligandsdf.is_file()):  
#         print(f"{proteinpdb} or {ligandsdf} does not exist.")  
#         return None  
      
#     # 生成口袋PDB文件（如果不存在）  
#     pocketpdb = proteinpdb.parent / (proteinpdb.name.rsplit('.', 1)[0] + '_pocket.pdb')  
#     if not pocketpdb.is_file():  
#         pymol_pocket(proteinpdb, ligandsdf, pocketpdb)  
      
#     # 跳过SDF转换步骤，直接使用PDB文件  
  
#     try:  
#         ligand = read_mol(ligandsdf)  
#         # 直接从PDB文件读取口袋信息  
#         pocket_dict = mol2graph_protein_from_pdb(pocketpdb)  
          
#         proinfo, liginfo = get_info(proteinpdb, ligandsdf)  
          
#         # 手动构建res字典  
#         ligand_dict = mol2graph_ligand(ligand)  
          
#         res = {  
#             'lc': ligand_dict['coords'], 'lf': ligand_dict['node_feat'],   
#             'lei': ligand_dict['edge_index'], 'lea': ligand_dict['edge_feat'],  
#             'pc': pocket_dict['coords'], 'pf': pocket_dict['node_feat'],   
#             'pei': pocket_dict['edge_index'], 'pea': pocket_dict['edge_feat'],  
#             'pdbid': name,  
#             'ligand_smiles': ligand_dict['smiles'],  
#             'protein_atom_names': pocket_dict['pro_name'],  
#             'protein_aa_names': pocket_dict['AA_name']  
#         }  
          
#         res['rfscore'] = RF_score(liginfo, proinfo)  
#         res['gbscore'] = GB_score(liginfo, proinfo)  
#         res['ecif'] = np.array(GetECIF(str(proteinpdb), str(ligandsdf)))  
          
#     except RuntimeError as e:  
#         print(proteinpdb, pocketpdb, ligandsdf, "Fail on reading molecule")  
#         return None  
  
#     # 其余处理逻辑保持不变  
#     ligand = (res['lc'], res['lf'], res['lei'], res['lea'])  
#     pocket = (res['pc'], res['pf'], res['pei'], res['pea'])  
      
#     try:  
#         raw = gen_graph(ligand, pocket, name, protein_cutoff=protein_cutoff, pocket_cutoff=pocket_cutoff, spatial_cutoff=spatial_cutoff)  
#     except ValueError as e:  
#         print(f"{name}: Error gen_graph from raw feature {str(e)}")  
#         return None  
      
#     # 返回字典格式  
#     result_dict = {  
#         'edge_index': raw[2], 'edge_feat': raw[3], 'node_feat': raw[1], 'coords': raw[0],  
#         'pro_name': res['protein_atom_names'], 'AA_name': res['protein_aa_names'],  
#         'smiles': res['ligand_smiles'], 'rmsd': rmsd,  
#         'rfscore': res['rfscore'], 'gbscore': res['gbscore'], 'ecif': res['ecif'],  
#         'pk': pk, 'pdbid': name, 'num_node': raw[4], 'num_edge': raw[5]  
#     }  
      
#     return result_dict

def parallel_helper(proteinpdb, ligandsdf, name, pk, rmsd, protein_cutoff, pocket_cutoff, spatial_cutoff):    
    """    
    处理单个蛋白质-配体复合物的并行辅助函数    
    集成PLIP分析到特征生成流程    
    """    
    RDLogger.DisableLog('rdApp.*')  # 关闭 RDKit 日志输出    
        
    if not (proteinpdb.is_file() and ligandsdf.is_file()):    
        print(f"{proteinpdb} or {ligandsdf} does not exist.")    
        return None    
        
    # 生成口袋PDB文件（如果不存在）    
    pocketpdb = proteinpdb.parent / (proteinpdb.name.rsplit('.', 1)[0] + '_pocket.pdb')    
    if not pocketpdb.is_file():    
        pymol_pocket(proteinpdb, ligandsdf, pocketpdb)    
        
    try:    
        ligand = read_mol(ligandsdf)    
        # 直接从PDB文件读取口袋信息    
        pocket_dict = mol2graph_protein_from_pdb(pocketpdb)    
            
        proinfo, liginfo = get_info(proteinpdb, ligandsdf)    
            
        # 手动构建res字典    
        ligand_dict = mol2graph_ligand(ligand)    
            
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
            
        # 计算评分特征    
        res['gbscore'] = GB_score(liginfo, proinfo)    
            
        # 新增：PLIP相互作用分析    
        plip_interactions = analyze_plip_interactions(str(proteinpdb), str(ligandsdf))    
        if plip_interactions:    
            print(f"发现的结合位点数: {len(plip_interactions)}")    
        res['plip_interactions'] = plip_interactions    
            
    except RuntimeError as e:    
        print(proteinpdb, pocketpdb, ligandsdf, "Fail on reading molecule")    
        return None    
        
    # 提取配体与 pocket 特征    
    ligand = (res['lc'], res['lf'], res['lei'], res['lea'])    
    pocket = (res['pc'], res['pf'], res['pei'], res['pea'])    
        
    try:    
        # 修改：传递PLIP分析结果、文件路径和原子ID映射信息给gen_graph    
        raw = gen_graph(  
            ligand, pocket, name,  
            protein_cutoff=protein_cutoff,  
            pocket_cutoff=pocket_cutoff,  
            spatial_cutoff=spatial_cutoff,  
            protein_file=str(proteinpdb),  
            ligand_file=str(ligandsdf),  
            plip_interactions=res['plip_interactions'],  
            ligand_dict=ligand_dict,  # 新增：传递配体字典（包含原子ID信息）  
            pocket_dict=pocket_dict   # 新增：传递蛋白字典（包含原子ID信息）  
        )  
  
        comp_coord, comp_feat, comp_ei, comp_ea, comp_num_node, comp_num_edge, lig_sei, lig_sea, pro_sei, pro_sea = raw  
    except ValueError as e:    
        print(f"{name}: Error gen_graph from raw feature {str(e)}")    
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
        'pro_spatial_edge_index': pro_sei + len(comp_feat) // 2,  # 蛋白空间边下标偏移  
        'pro_spatial_edge_attr': pro_sea  
    }  
        
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