import argparse  # 命令行参数解析
from pathlib import Path  # 跨平台路径对象
from tqdm import tqdm  # 进度条（本脚本未用到）
import pickle  # 对象序列化/反序列化
from preprocess import gen_feature, gen_graph, to_pyg_graph, get_info,  GB_score,  analyze_plip_interactions  # 自定义化学特征/图处理函数
from joblib import Parallel, delayed  # 并行处理
from utils import read_mol, obabel_pdb2mol, pymol_pocket, correct_sanitize_v2  # 分子文件处理工具
import numpy as np  # 数值处理
from rdkit import Chem, RDLogger  # RDKit 化学库及日志
import tempfile  # 临时文件操作
import pandas as pd  # 表格处理（本脚本未用到）
import os  # 操作系统相关，如文件删除


# 多线程单个样本处理主函数
def parallel_helper(proteinpdb, ligandsdf, name_prefix, mol, pdict, protein_cutoff, pocket_cutoff, spatial_cutoff):  
    RDLogger.DisableLog('rdApp.*')  # 关闭 RDKit 日志输出，防止多线程干扰  
    # 创建临时文件用于中间数据保存  
    _, templigand = tempfile.mkstemp(suffix='.sdf')  # 临时配体 sdf  
    os.close(_)  # 关闭文件描述符  
    _, temppocketpdb = tempfile.mkstemp(suffix='.pdb')  # 临时 pocket pdb  
    os.close(_)  
    _, temppocketsdf = tempfile.mkstemp(suffix='.sdf')  # 临时 pocket sdf  
    os.close(_)  
      
    pymol_pocket(proteinpdb, ligandsdf, temppocketpdb)  # 用 pymol 选出 pocket 区域保存为 pdb  
    obabel_pdb2mol(temppocketpdb, temppocketsdf)  # 用 openbabel 转为 sdf 格式  
    assert "_Name" in pdict, f'Property dict should have _Name key, but currently: {pdict}'  # 检查分子属性有名字  
    name = name_prefix + f'_{pdict["_Name"]}'  # 生成唯一名字  
      
    try:  
        ligand = correct_sanitize_v2(mol)  # RDKit 分子消毒修正  
        Chem.MolToMolFile(ligand, templigand)  # 保存修正结果为 sdf  
        pocket = read_mol(temppocketsdf)  # 读取 pocket 分子  
        proinfo, liginfo = get_info(proteinpdb, templigand)  # 获取蛋白与配体信息（如原子坐标、类型等）  
        res = gen_feature(ligand, pocket, name)  # 生成结构特征（如原子、键、环境等）  
          
        # 机器学习相关评分特征  
        # res['rfscore'] = RF_score(liginfo, proinfo)  # 随机森林打分  
        res['gbscore'] = GB_score(liginfo, proinfo)  # 梯度提升树打分   
          
        # 新增：PLIP相互作用分析  
        print(f"正在进行PLIP分析: {name}")  
        plip_interactions = analyze_plip_interactions(str(proteinpdb), str(templigand))  
        res['plip_interactions'] = plip_interactions  
          
    except RuntimeError as e:  
        # 只要分子读入失败就报错并返回 None  
        print(proteinpdb, temppocketsdf, templigand, "Fail on reading molecule")  
        return None  
  
    # 提取配体与 pocket 特征（如原子类别、坐标、环境等，作为图节点/边属性）  
    ligand = (res['lc'], res['lf'], res['lei'], res['lea'])  
    pocket = (res['pc'], res['pf'], res['pei'], res['pea'])  
  
    # 生成分子复合物图，异常时跳过  
    try:  
        # 修改：传递PLIP分析结果和文件路径给gen_graph  
        raw = gen_graph(  
            ligand, pocket, name,   
            protein_cutoff=protein_cutoff,  
            pocket_cutoff=pocket_cutoff,  
            spatial_cutoff=spatial_cutoff,  
            protein_file=str(proteinpdb),  
            ligand_file=str(templigand),  
            plip_interactions=res['plip_interactions']  
        )  # 生成原始图结构特征（如邻接、节点、边等），cutoff 控制边界  
    except ValueError as e:  
        print(f"{name}: Error gen_graph from raw feature {str(e)}")  
        return None  
  
    # 转换为 PyTorch Geometric 图数据结构，并加上全局特征  
    graph = to_pyg_graph(  
        list(raw) + [ res['gbscore'], res['ecif'], -1, name],   
        frame=-1, rmsd_lig=0.0, rmsd_pro=0.0  
    )  
  
    # 清理所有临时文件，防止磁盘爆满  
    os.remove(templigand)  
    os.remove(temppocketpdb)  
    os.remove(temppocketsdf)  
    return graph  # 返回该分子复合物的图对象

# 并行批量处理复合物
def process_complex(proteinpdb: Path, ligandsdf: Path, name_prefix: str, njobs: int, protein_cutoff, pocket_cutoff, spatial_cutoff):
    suppl = Chem.SDMolSupplier(str(ligandsdf), sanitize=False, strictParsing=False)  # 读取 ligand sdf 文件，允许多 conformer
    mols = list(suppl)  # 转成分子列表
    graphs = []

    # joblib 并行分发任务（每个 ligand 一个进程）
    res = Parallel(n_jobs=njobs)(
        delayed(parallel_helper)(
            proteinpdb, ligandsdf, f"{idx}_{name_prefix}", mol, mol.GetPropsAsDict(True),
            protein_cutoff, pocket_cutoff, spatial_cutoff
        )
        for idx, mol in enumerate(mols)
    )
    # 收集所有图对象
    for i in res:
        if i: graphs.append(i)
    return graphs  # 返回所有图对象列表


if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # 命令行参数解析器
    parser.add_argument('protein', type=Path)  # 蛋白 pdb 文件
    parser.add_argument('ligand', type=Path)  # 配体 sdf/mol 文件
    parser.add_argument('output', type=Path)  # 输出文件名（pkl 格式保存所有图的列表）
    parser.add_argument('--njobs', type=int, default=8)  # 并发进程数，-1 用全部 CPU
    parser.add_argument('--protein_cutoff', type=float, default=5.)  # 蛋白原子边界距离（Å）
    parser.add_argument('--pocket_cutoff', type=float, default=5.)  # pocket 原子边界距离（Å）
    parser.add_argument('--spatial_cutoff', type=float, default=5.)  # 图生成时空间距离阈值

    args = parser.parse_args()  # 解析参数

    # 文件格式检查
    if args.protein.name.split('.')[-1] != 'pdb':
        raise ValueError('Make sure your protein file is in pdb format.')
    if args.ligand.name.split('.')[-1] not in ['sdf', 'mol']:
        raise ValueError("Make sure your ligand file is in sdf/mol format.")

    protein, ligand, output = args.protein, args.ligand, args.output

    name_prefix = protein.name.rsplit('.', 1)[0]  # 去掉后缀作为样本名前缀

    # 批量处理，得到所有复合物的图
    graphs = process_complex(
        protein, ligand, name_prefix, args.njobs,
        args.protein_cutoff, args.pocket_cutoff, args.spatial_cutoff
    )

    # 所有图对象序列化保存为 pickle 文件
    with open(args.output, 'wb') as f:
        pickle.dump(graphs, f)