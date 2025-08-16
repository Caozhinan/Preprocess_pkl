import argparse
from multiprocessing.sharedctypes import Value
from pathlib import Path
from tqdm import tqdm, trange
from sklearn.model_selection import train_test_split
import pickle
from preprocess import gen_feature, load_pk_data, gen_graph, to_pyg_graph, get_info,  GB_score
from joblib import Parallel, delayed
from utils import read_mol
import pandas as pd
import numpy as np

# 并行处理每个分子动力学采样条目
def parallel_helper(refined_path: Path, name: str):
    frames = []
    # 读取当前样本的 RMSD 数据
    rmsd_df = pd.read_csv(refined_path / name / "rmsd.csv")
    for i in range(100):  # 采样 100 帧
        try:
            # 读取配体和口袋的 SDF 文件
            ligand = read_mol(refined_path / name / f"{name}_ligand_{i}.sdf")
            pocket = read_mol(refined_path / name / f"{name}_pocket_{i}.sdf")
            # 生成配体-口袋特征
            res = gen_feature(ligand, pocket, name)
            res['frame'] = i
            # 读取对应帧的 RMSD 信息
            res['rmsd_lig'] = rmsd_df[rmsd_df['frame'] == i]['rmsd_lig'].item()
            res['rmsd_pro'] = rmsd_df[rmsd_df['frame'] == i]['rmsd_pro'].item()
            # 获取蛋白和配体的详细信息
            proinfo, liginfo = get_info(
                refined_path / name / f"{name}_protein_{i}.pdb",
                refined_path / name / f"{name}_ligand_{i}.sdf"
            )
            # 计算各种打分
            # res['rfscore'] = RF_score(liginfo, proinfo)
            res['gbscore'] = GB_score(liginfo, proinfo)
        except RuntimeError as e:
            continue  # 跳过异常帧
        except OSError as e:
            continue  # 跳过异常帧
        frames.append(res)
    return frames  # 返回该条目的所有帧的特征

# 主数据处理函数
def process_dynamics(
    core_path: Path,
    refined_path: Path,
    md_path: Path,
    dataset_name: str,
    output_path: Path,
    protein_cutoff: float,
    pocket_cutoff: float,
    spatial_cutoff: float,
    seed: int
):
    # 读取核心集、精炼集、MD采样集的所有条目名
    core_set_list = [d.name for d in core_path.iterdir() if len(d.name) == 4]
    refined_set_list = [d.name for d in refined_path.iterdir() if len(d.name) == 4]
    md_set_list = [d.name for d in md_path.iterdir()]
    
    # 读取结合常数数据（pK），合并为字典
    pk_file_gen = list((refined_path / 'index').glob('INDEX_general_PL_data.*'))[0]
    pk_file_ref = list((refined_path / 'index').glob('INDEX_refined_data.*'))[0]
    print(f"Loading pk data from {pk_file_gen} and {pk_file_ref}")
    pk_dict = {**load_pk_data(pk_file_gen), **load_pk_data(pk_file_ref)}
    # 检查所有核心集、精炼集的PDB ID都在pk_dict里
    assert set(pk_dict.keys()).issuperset(refined_set_list)
    assert set(pk_dict.keys()).issuperset(core_set_list)
    print(f"Total {len(md_set_list)} md data")
    # 特征生成与缓存
    if not Path('./rawgraph.pkl').is_file():
        # 并行处理所有 md_set_list（每个分子模型，100帧），22核加速
        res = Parallel(n_jobs=22)(
            delayed(parallel_helper)(md_path, name)
            for name in tqdm(md_set_list, desc="Load refined", ncols=80)
        )
        pickle.dump(res, open('./rawgraph.pkl', 'wb'))  # 缓存结果
    else:
        res = pickle.load(open('./rawgraph.pkl', 'rb'))  # 直接加载缓存

    # 按 pdbid 组织特征
    processed_dict = {feat[0]['pdbid']: feat for feat in res}
    # 补充 pK 数据
    for pdbid in processed_dict:
        if pdbid[:4] not in pk_dict:
            print("pdbid:", pdbid, "not in pk_dict")
        pk = pk_dict[pdbid[:4]]
        for i in range(len(processed_dict[pdbid])):
            processed_dict[pdbid][i]['pk'] = pk

    refined_data, core_data, logstr = [], [], []
    # 遍历每个条目（pdbid）及其所有帧
    for name, vv in tqdm(list(processed_dict.items()), desc="Construct graph", ncols=80):
        graphs = []
        for v in tqdm(vv, ncols=80, leave=False, desc="Frame"):
            # 取出配体和口袋的特征，用于生成图
            ligand = (v['lc'], v['lf'], v['lei'], v['lea'])
            pocket = (v['pc'], v['pf'], v['pei'], v['pea'])
            try:
                # 生成蛋白-配体图结构，返回边数量等
                raw = gen_graph(
                    ligand, pocket, name,
                    protein_cutoff=protein_cutoff,
                    pocket_cutoff=pocket_cutoff,
                    spatial_cutoff=spatial_cutoff
                )
                edge_nums = raw[5]
                # 筛选：蛋白节点数、蛋白-配体连接数太少的丢弃
                if edge_nums[1] <= 3:
                    raise ValueError("<4 protein edges (fewer nodes)")
                if edge_nums[2] <= 3:
                    raise ValueError("<4 protein-ligand edges (fewer nodes)")
            except ValueError as e:
                logstr.append(f"{name} - {v['frame']}: Error gen_graph from raw feature: {str(e)}")
                continue
            except IndexError as e:
                logstr.append(f"{name} - {v['frame']}: Error gen_graph from raw feature: {str(e)}")
                continue
            # 转为 pyg 格式图（用于深度学习）
            g = to_pyg_graph(
                list(raw) + [v['rfscore'], v['gbscore'], v['ecif'], v['pk'], name],
                frame=v['frame'], rmsd_lig=v['rmsd_lig'], rmsd_pro=v['rmsd_pro']
            )
            # 距离筛查，最大原子距离超过 45A 的丢弃
            m = (g['pos'].unsqueeze(0) - g['pos'].unsqueeze(1)).norm(dim=-1).max()
            if m > 45:
                logstr.append(f"{name} - {v['frame']}: max dist {m:.3f}, maybe out of boundary")
                continue
            graphs.append(g)
        # 按是否在核心集分为 test/train
        if name[:4] in core_set_list:
            core_data.append(graphs)
        else:
            refined_data.append(graphs)

    print(f"Got {len(refined_data)} / {len(refined_set_list)} graphs for train and val, {len(core_data)} / {len(core_set_list)} for test")
    # 输出到指定目录
    if not output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)
    with open(output_path / f"{dataset_name}_train_val.pkl", 'wb') as f:
        pickle.dump(refined_data, f)
    with open(output_path / f"{dataset_name}_test.pkl", 'wb') as f:
        pickle.dump(core_data, f)
    print(*logstr, sep='\n')  # 打印日志

if __name__ == "__main__":
    # 命令行参数解析
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path_core', type=Path)
    parser.add_argument('data_path_refined', type=Path)
    parser.add_argument('md_path', type=Path)
    parser.add_argument('output_path', type=Path)
    parser.add_argument('--dataset_name', type=str, default='pdbbind')
    parser.add_argument('--protein_cutoff', type=float, default=5.)
    parser.add_argument('--pocket_cutoff', type=float, default=5.)
    parser.add_argument('--spatial_cutoff', type=float, default=5.)
    parser.add_argument('--seed', type=int, default=2022)

    args = parser.parse_args()
    process_dynamics(
        args.data_path_core, args.data_path_refined, args.md_path,
        args.dataset_name, args.output_path,
        args.protein_cutoff, args.pocket_cutoff, args.spatial_cutoff, args.seed
    )

# 2. parallel_helper 函数

#     输入：数据存储路径和条目名（如 PDB ID）。
#     遍历指定路径下的 100 个（0~99）分子动力学采样帧，对每一帧：
#         读取配体和口袋的结构文件（SDF）。
#         生成特征（如分子描述符）。
#         读取并关联 RMSD 信息。
#         计算打分（RFScore、GBScore、ECIF）。
#         错误跳过，成功则加入结果列表。
#     返回所有帧的特征列表。

# 3. process_dynamics 主处理函数

#     读取核心集、精炼集和MD采样集的所有条目（目录名）。
#     读取结合常数（pK）数据，合并为字典。
#     检查所有精炼集、核心集的 PDB ID 都在 pK 数据里。
#     如果没有原始特征数据缓存（rawgraph.pkl），则多进程并行计算所有 MD 采样集的特征，保存为 pkl；否则直接加载。
#     按 pdbid 组织特征，将对应的 pK 数据填充进去。
#     遍历每个条目的所有帧，尝试生成图结构（pyg），并做合理性检查（如边数、最大原子间距离）。
#     按照核心集和精炼集将数据分为 test 和 train/val 集。
#     将处理好的数据使用 pickle 序列化保存到输出目录。

# 4. 主程序入口

#     argparse 解析命令行参数，调用主处理函数。
