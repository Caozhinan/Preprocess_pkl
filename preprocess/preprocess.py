
import numpy as np  
from scipy.spatial.distance import cdist  
from pathlib import Path  
# from torch_geometric.data import Data  
# import torch  
from rdkit import Chem  
from pymol import cmd  
# 添加PLIP导入  
from plip.structure.preparation import PDBComplex  
import tempfile  
import os  
from plip_analysis import merge_protein_ligand_with_pymol  
from intra_pro_plip import (  
    InteractionAnalyzer,  
    AtomInfo, Config )
from water_metal_detection import detect_water_bridges_from_atoms, detect_metal_complex_from_atoms
from intra_lig_plip import InteractionAnalyzer as LigandInteractionAnalyzer
import time
from openbabel import pybel

# 扩展边特征编码  
SPATIAL_EDGE = [4, 0, 0]  # 原有空间边  
HYDROGEN_BOND_EDGE = [5, 1, 0]  # 氢键  
HYDROPHOBIC_EDGE = [5, 2, 0]    # 疏水相互作用  
PI_STACKING_EDGE = [5, 3, 0]    # π-π堆积  
PI_CATION_EDGE = [5, 4, 0]      # π-阳离子相互作用  
SALT_BRIDGE_EDGE = [5, 5, 0]    # 盐桥  
WATER_BRIDGE_EDGE = [5, 6, 0]   # 水桥  
HALOGEN_BOND_EDGE = [5, 7, 0]   # 卤键  
METAL_COMPLEX_EDGE = [5, 8, 0]  # 金属配位  
OTHERS_EDGE = [5, 9, 0]         # 其他相互作用  
  
# 相互作用类型映射  
INTERACTION_TYPE_MAP = {  
    'hydrogen_bonds': HYDROGEN_BOND_EDGE,        # 改为复数  
    'hydrophobic_contacts': HYDROPHOBIC_EDGE,    # 改为复数  
    'pi_stacking': PI_STACKING_EDGE,  
    'pi_cation': PI_CATION_EDGE,  
    'salt_bridges': SALT_BRIDGE_EDGE,            # 改为复数  
    'water_bridges': WATER_BRIDGE_EDGE,          # 改为复数  
    'halogen_bonds': HALOGEN_BOND_EDGE,          # 改为复数  
    'metal_complexes': METAL_COMPLEX_EDGE,       # 改为复数  
    'others': OTHERS_EDGE  
}

def analyze_plip_interactions(protein_file, ligand_file):  
    """  
    使用PLIP分析蛋白质-配体相互作用  
      
    Args:  
        protein_file: 蛋白质PDB文件路径  
        ligand_file: 配体SDF文件路径  
      
    Returns:  
        dict: 包含原子对相互作用信息的字典  
    """  
    # Create complex file in the same directory as the protein file  
    protein_path = Path(protein_file)  
    complex_file = protein_path.parent / "complex.pdb"  
      
     
    if not complex_file.exists():  
        # print(f"Creating complex file: {complex_file}")  
          
        if not merge_protein_ligand_with_pymol(protein_file, ligand_file, str(complex_file)):  
            print(f"Failed to create complex file")  
            return {}  
          
        # 新增：修改HETATM行的链标识符  
        fix_hetatm_chain_ids(str(complex_file), 'X')  
        
    try:  
        my_mol = PDBComplex()  
        my_mol.load_pdb(str(complex_file))
        # print(f"DEBUG: 识别出的配体数量: {len(my_mol.ligands)}")  
        # for i, ligand in enumerate(my_mol.ligands):  
        #     print(f"DEBUG: 配体 {i+1}: {ligand.hetid}, 链: {ligand.chain}, 位置: {ligand.position}")  
            # print(f"DEBUG: 配体类型: {ligand.type}, 长名称: {ligand.longname}")  

        # print(f"DEBUG: 排除的分子: {my_mol.excluded}") 
        my_mol.analyze()  
          
        # 提取原子对信息  
        atom_pairs_dict = {}  
        # print(f"DEBUG: 发现的结合位点标识符: {list(my_mol.interaction_sets.keys())}")  
        # print(f"DEBUG: 结合位点详细信息:")
        for bsid, interactions in my_mol.interaction_sets.items():  
            atom_pairs_dict[bsid] = {  
                'hydrogen_bonds': [],  
                'hydrophobic_contacts': [],  
                'pi_stacking': [],  
                'pi_cation': [],  
                'salt_bridges': [],  
                'water_bridges': [],  
                'halogen_bonds': [],  
                'metal_complexes': []  
            }  
              
            # 氢键  
            for hbond in interactions.hbonds_ldon + interactions.hbonds_pdon:  
                atom_pairs_dict[bsid]['hydrogen_bonds'].append({  
                    'protein_coords': hbond.d.coords if hbond.protisdon else hbond.a.coords,  
                    'ligand_coords': hbond.a.coords if hbond.protisdon else hbond.d.coords,  
                    'distance': hbond.distance_ad  
                })  
              
            # 疏水相互作用  
            for hydrophobic in interactions.hydrophobic_contacts:  
                atom_pairs_dict[bsid]['hydrophobic_contacts'].append({  
                    'protein_coords': hydrophobic.bsatom.coords,  
                    'ligand_coords': hydrophobic.ligatom.coords,  
                    'distance': hydrophobic.distance  
                })  
              
            # π-π堆积  
            for pistack in interactions.pistacking:  
                atom_pairs_dict[bsid]['pi_stacking'].append({  
                    'protein_coords': pistack.proteinring.center,  
                    'ligand_coords': pistack.ligandring.center,  
                    'distance': pistack.distance  
                })  
              
            # π-阳离子相互作用  
            for pication in interactions.pication_laro + interactions.pication_paro:  
                atom_pairs_dict[bsid]['pi_cation'].append({  
                    'protein_coords': pication.ring.center if pication.protcharged else pication.charge.center,  
                    'ligand_coords': pication.charge.center if pication.protcharged else pication.ring.center,  
                    'distance': pication.distance  
                })  
              
            # 盐桥  
            for saltbridge in interactions.saltbridge_lneg + interactions.saltbridge_pneg:  
                atom_pairs_dict[bsid]['salt_bridges'].append({  
                    'protein_coords': saltbridge.positive.center if saltbridge.protispos else saltbridge.negative.center,  
                    'ligand_coords': saltbridge.negative.center if saltbridge.protispos else saltbridge.positive.center,  
                    'distance': saltbridge.distance  
                })  
              
            # 水桥  
            for wbridge in interactions.water_bridges:  
                atom_pairs_dict[bsid]['water_bridges'].append({  
                    'protein_coords': wbridge.d.coords if wbridge.protisdon else wbridge.a.coords,  
                    'ligand_coords': wbridge.a.coords if wbridge.protisdon else wbridge.d.coords,  
                    'distance': (wbridge.distance_dw + wbridge.distance_aw) / 2  # 平均距离  
                })  
              
            # 卤键  
            for halogen in interactions.halogen_bonds:  
                atom_pairs_dict[bsid]['halogen_bonds'].append({  
                    'protein_coords': halogen.acc.coords,  
                    'ligand_coords': halogen.don.coords,  
                    'distance': halogen.distance  
                })  
              
            # 金属配位  
            for metal in interactions.metal_complexes:  
                atom_pairs_dict[bsid]['metal_complexes'].append({  
                    'protein_coords': metal.target.atom.coords,  
                    'ligand_coords': metal.metal.coords,  
                    'distance': metal.distance  
                })  
          
        return atom_pairs_dict  
          
    except Exception as e:  
        print(f"PLIP分析失败: {e}")  
        return {}  

def fix_hetatm_chain_ids(pdb_file, chain_id='X'):  
    """  
    修改PDB文件中HETATM行的链标识符  
      
    Args:  
        pdb_file: PDB文件路径  
        chain_id: 要设置的链标识符，默认为'X'  
    """  
    import os  
      
    # 读取原文件  
    with open(pdb_file, 'r') as f:  
        lines = f.readlines()  
      
    # 修改HETATM行  
    modified_lines = []  
    for line in lines:  
        if line.startswith('HETATM'):  
            # PDB格式：HETATM行的链标识符在第22位（索引21）  
            # 原格式：HETATM   48  H   UNK     0     105.436...  
            # 目标格式：HETATM   48  H   UNK X   0     105.436...  
            if len(line) > 21:  
                # 将第22位（索引21）设置为链标识符  
                line = line[:21] + chain_id + line[22:]  
            modified_lines.append(line)  
        else:  
            modified_lines.append(line)  
      
    # 写回文件  
    with open(pdb_file, 'w') as f:  
        f.writelines(modified_lines)  
      
    # print(f"✓ 已修改PDB文件中的HETATM链标识符为: {chain_id}")

def convert_plip_to_edges(atom_pairs_dict, lig_coord, pro_coord, lig_num_atom,   
                         lig_atom_id_map=None, pro_atom_id_map=None,   
                         lig_mol=None, pro_mol=None):  
    """  
    将PLIP相互作用转换为图边  
      
    Args:  
        atom_pairs_dict: PLIP分析结果  
        lig_coord: 配体原子坐标  
        pro_coord: 蛋白原子坐标  
        lig_num_atom: 配体原子数量  
        lig_atom_id_map: 配体原子ID到图索引的映射 (新增)  
        pro_atom_id_map: 蛋白原子ID到图索引的映射 (新增)  
        lig_mol: 配体分子对象 (新增)  
        pro_mol: 蛋白分子对象 (新增)  
      
    Returns:  
        tuple: (edge_index, edge_attr) 包含相互作用类型和距离的边  
    """  
    edge_list = []  
    edge_attr_list = []  
      
    for bsid, interactions in atom_pairs_dict.items():  
        for interaction_type, pairs in interactions.items():  
            if not pairs:  
                continue  
                  
            edge_type = INTERACTION_TYPE_MAP.get(interaction_type, OTHERS_EDGE)  
              
            for pair in pairs:  
                # 优先使用原子ID进行查找  
                lig_idx = None  
                pro_idx = None  
                  
                # 如果有原子ID信息，直接使用ID查找  
                if (lig_atom_id_map is not None and pro_atom_id_map is not None and   
                    'ligand_atom_id' in pair and 'protein_atom_id' in pair):  
                      
                    lig_atom_id = pair['ligand_atom_id']  
                    pro_atom_id = pair['protein_atom_id']  
                      
                    if lig_atom_id in lig_atom_id_map:  
                        lig_idx = lig_atom_id_map[lig_atom_id]  
                    if pro_atom_id in pro_atom_id_map:  
                        pro_idx = pro_atom_id_map[pro_atom_id]  
                  
                # 如果ID查找失败，回退到坐标查找  
                if lig_idx is None or pro_idx is None:  
                    lig_distances = cdist([pair['ligand_coords']], lig_coord)[0]  
                    pro_distances = cdist([pair['protein_coords']], pro_coord)[0]  
                      
                    lig_idx = np.argmin(lig_distances)  
                    pro_idx = np.argmin(pro_distances)  
                  
                # 实际欧氏距离  
                actual_distance = np.linalg.norm(  
                    np.array(pair['ligand_coords']) - np.array(pair['protein_coords'])  
                )  
                  
                # 创建边特征：[类型编码] + [欧氏距离]  
                edge_feature = edge_type + [actual_distance]  
                  
                # 添加双向边  
                edge_list.append([lig_idx, pro_idx + lig_num_atom])  
                edge_list.append([pro_idx + lig_num_atom, lig_idx])  
                edge_attr_list.append(edge_feature)  
                edge_attr_list.append(edge_feature)  
      
    if edge_list:  
        edge_index = np.array(edge_list, dtype=np.int64).T  
        edge_attr = np.array(edge_attr_list, dtype=np.float32)  
    else:  
        edge_index = np.empty((2, 0), dtype=np.int64)  
        edge_attr = np.empty((0, 4), dtype=np.float32)  
      
    return edge_index, edge_attr

def dim2(arr: np.ndarray):
    # 判断数组是否为二维
    return len(arr.shape) == 2

def reindex(atom_idx: list, edge_index: np.ndarray):
    # 对边的节点编号重新映射（因为去除部分原子后索引变化）
    if len(edge_index.shape) != 2:
        return edge_index
    indexmap = {old: new for new, old in enumerate(atom_idx)}
    mapfunc = np.vectorize(indexmap.get)
    edge_index_new = np.array([mapfunc(edge_index[0]), mapfunc(edge_index[1])])
    return edge_index_new

def remove_duplicated_edges(ei: np.ndarray, ea: np.ndarray, ref_ei: np.ndarray):
    # 删除和已有边重复的边（避免冗余）
    if not len(ea):
        return ei, ea
    existing_set = {(i, j) for i, j in zip(*ref_ei)}
    mask = []
    for i, j in zip(*ei):
        if (i, j) in existing_set:
            mask.append(False)  # 已有则删除
        else:
            mask.append(True)  # 否则保留
    mask = np.array(mask, dtype=bool)
    ei_n = np.array([ei[0][mask], ei[1][mask]])
    ea_n = ea[mask]
    return ei_n, ea_n

# def gen_feature(ligand: Chem.Mol, pocket: Chem.Mol, pdbid: str):  
#     ligand_dict = mol2graph_ligand(ligand)  # 使用小分子版本  
#     pocket_dict = mol2graph_protein_from_pdb(pocket)  # 使用蛋白质版本  
      
#     ligand_coords, ligand_features, ligand_edge_index, ligand_edge_attr = ligand_dict['coords'], ligand_dict[  
#         'node_feat'], ligand_dict['edge_index'], ligand_dict['edge_feat']  
#     pocket_coords, pocket_features, pocket_edge_index, pocket_edge_attr = pocket_dict['coords'], pocket_dict[  
#         'node_feat'], pocket_dict['edge_index'], pocket_dict['edge_feat']  
  
#     # 形状检查  
#     if not (dim2(ligand_coords) and dim2(ligand_features) and dim2(ligand_edge_index) and dim2(ligand_edge_attr)):  
#         raise RuntimeError(f"Ligand feature shape error")  
#     if not (dim2(pocket_coords) and dim2(pocket_features) and dim2(pocket_edge_index) and dim2(pocket_edge_attr)):  
#         raise RuntimeError(f"Protein feature shape error")  
  
#     return {'lc': ligand_coords, 'lf': ligand_features, 'lei': ligand_edge_index, 'lea': ligand_edge_attr,  
#             'pc': pocket_coords, 'pf': pocket_features, 'pei': pocket_edge_index, 'pea': pocket_edge_attr,  
#             'pdbid': pdbid,  
#             'ligand_smiles': ligand_dict['smiles'],  # 新增  
#             'protein_atom_names': pocket_dict['pro_name'],  # 新增  
#             'protein_aa_names': pocket_dict['AA_name']}  # 新增

def gen_spatial_edge(dm: np.ndarray, spatial_cutoff: float = 5):  
    # 生成所有空间距离<cutoff的节点对，并赋空间边特征  
    if spatial_cutoff <= 0.1:  
        return np.empty((2, 0), dtype=np.int64), np.empty((0, 4), dtype=np.float32)  
      
    src, dst = np.where((dm <= spatial_cutoff) & (dm > 0.1))  
    edge_index = [(x, y) for x, y in zip(src, dst)]  
      
    # 为每条边添加距离特征  
    edge_attr = []  
    for x, y in edge_index:  
        distance = dm[x, y]  
        edge_feature = SPATIAL_EDGE + [distance]  # [4, 0, 0, distance]  
        edge_attr.append(edge_feature)  
  
    edge_index = np.array(edge_index, dtype=np.int64).T  
    edge_attr = np.array(edge_attr, dtype=np.float32)  # 改为float32  
    return edge_index, edge_attr

def gen_ligpro_edge(dm: np.ndarray, pocket_cutoff: float, plip_interactions=None,   
                   lig_coord=None, pro_coord=None, lig_atom_id_map=None, pro_atom_id_map=None):  
    # print(f"DEBUG: plip_interactions = {plip_interactions}")    
    # print(f"DEBUG: plip_interactions is None = {plip_interactions is None}")    
        
    # if plip_interactions:    
    #     print(f"DEBUG: Found {len(plip_interactions)} binding sites")    
    #     for bsid, interactions in plip_interactions.items():    
    #         print(f"DEBUG: Binding site {bsid} has interactions: {list(interactions.keys())}")    
      
    # 生成配体-蛋白之间的边，包括PLIP相互作用和空间边    
    lig_num_atom, pro_num_atom = dm.shape    
        
    # 如果有PLIP分析结果，使用convert_plip_to_edges函数    
    if plip_interactions and lig_coord is not None and pro_coord is not None:    
        return convert_plip_to_edges(  
            plip_interactions, lig_coord, pro_coord, lig_num_atom,  
            lig_atom_id_map=lig_atom_id_map,   
            pro_atom_id_map=pro_atom_id_map  
        )    
        
    # 否则生成普通空间边，但要包含距离特征    
    lig_idx, pro_idx = np.where(dm <= pocket_cutoff)    
    edge_list = []    
    edge_attr_list = []    
        
    for x, y in zip(lig_idx, pro_idx):    
        distance = dm[x, y]    
        edge_feature = SPATIAL_EDGE + [distance]  # [4, 0, 0, distance]    
            
        # 双向边    
        edge_list.extend([(x, y + lig_num_atom), (y + lig_num_atom, x)])    
        edge_attr_list.extend([edge_feature, edge_feature])    
        
    edge_index = np.array(edge_list, dtype=np.int64).T if edge_list else np.empty((2, 0), dtype=np.int64)    
    edge_attr = np.array(edge_attr_list, dtype=np.float32) if edge_attr_list else np.empty((0, 4), dtype=np.float32)    
        
    return edge_index, edge_attr
 
  
def classify_protein_spatial_edges(pro_sei, pro_sea, pro_coord, protein_file,     
                                 pro_atom_ids=None, pro_atom_id_map=None):  
    import time  
    if len(pro_sei) == 0:  
        return pro_sei, pro_sea  
      
    try:    
        molecule = pybel.readfile("pdb", protein_file).__next__()  
        analyzer = InteractionAnalyzer()  
        analyzer.set_molecule(molecule)  
          
        # 方案1：预构建坐标到原子的哈希映射（高效且准确）  
        coord_to_atom_map = {}  
        tolerance = 0.01  # 0.01Å精度，更精确的匹配  
          
        for atom in molecule.atoms:  
            coord_key = tuple(np.round(atom.coords, 2))  # 保留2位小数精度  
            coord_to_atom_map[coord_key] = atom  
          
        def find_atom_by_coord_fast(coord):  
            coord_key = tuple(np.round(coord, 2))  
            return coord_to_atom_map.get(coord_key)  
          
        classified_edge_attr = []  
        successful_matches = 0  
        failed_matches = 0  
        interaction_detected = 0  
        total_edges = len(pro_sei.T)  
          
        t_big_loop = time.perf_counter()  
        for i, (src_idx, tgt_idx) in enumerate(pro_sei.T):  
            src_coord = pro_coord[src_idx]  
            tgt_coord = pro_coord[tgt_idx]  
            distance = np.linalg.norm(src_coord - tgt_coord)  
              
            # 使用快速坐标查找  
            src_ob_atom = find_atom_by_coord_fast(src_coord)  
            tgt_ob_atom = find_atom_by_coord_fast(tgt_coord)  
              
            src_atom = AtomInfo(src_ob_atom) if src_ob_atom else None  
            tgt_atom = AtomInfo(tgt_ob_atom) if tgt_ob_atom else None  
              
            if src_atom and tgt_atom:  
                successful_matches += 1  
                interaction_result = analyzer.analyze_atom_pair(src_atom, tgt_atom, distance)  
                  
                if interaction_result and interaction_result['interaction_types']:  
                    interaction_detected += 1  
                    interaction_type = interaction_result['interaction_types'][0]  
                      
                    # if interaction_detected <= 10:  
                    #     print(f"[DEBUG] 边 {i}: 检测到相互作用 {interaction_type}, 距离: {distance:.2f}Å")  
                      
                    if interaction_type == 'hydrogen_bond':  
                        edge_type_name = 'hydrogen_bonds'  
                    elif interaction_type == 'hydrophobic':  
                        edge_type_name = 'hydrophobic_contacts'  
                    elif interaction_type == 'salt_bridge':  
                        edge_type_name = 'salt_bridges'  
                    elif interaction_type == 'halogen_bond':  
                        edge_type_name = 'halogen_bonds'  
                    else:  
                        edge_type_name = 'others'  
                else:  
                    edge_type_name = check_ring_interactions_improved(  
                        src_atom, tgt_atom, distance, analyzer.rings, analyzer.charged_groups  
                    )  
                    if edge_type_name == 'others' and detect_water_bridges_from_atoms(src_atom, tgt_atom, distance, analyzer.molecule):  
                        edge_type_name = 'water_bridges'  
                    if edge_type_name == 'others' and detect_metal_complex_from_atoms(src_atom, tgt_atom, distance, analyzer.molecule):  
                        edge_type_name = 'metal_complexes'  
            else:  
                failed_matches += 1  
                edge_type_name = 'others'  
              
            edge_type = INTERACTION_TYPE_MAP.get(edge_type_name, OTHERS_EDGE)  
            edge_feature = edge_type + [distance]  
            classified_edge_attr.append(edge_feature)  
          
        # print(f"[DEBUG] 原子匹配统计:")  
        # print(f"  成功匹配: {successful_matches}/{total_edges} ({successful_matches/total_edges*100:.1f}%)")  
        # print(f"  匹配失败: {failed_matches}/{total_edges} ({failed_matches/total_edges*100:.1f}%)")  
        # print(f"  检测到特定相互作用: {interaction_detected}/{successful_matches} ({interaction_detected/successful_matches*100:.1f}% of successful matches)")  
          
        # print(f"[classify_protein_spatial_edges] total big loop time: {time.perf_counter() - t_big_loop:.3f}s")  
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)  
          
    except Exception as e:  
        print(f"使用 intra_pro_plip 分析失败: {e}")  
        classified_edge_attr = []  
        for i, (src_idx, tgt_idx) in enumerate(pro_sei.T):  
            distance = np.linalg.norm(pro_coord[src_idx] - pro_coord[tgt_idx])  
            edge_feature = OTHERS_EDGE + [distance]  
            classified_edge_attr.append(edge_feature)  
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)  
      
    return pro_sei, classified_edge_attr

def classify_ligand_spatial_edges(lig_sei, lig_sea, lig_coord, ligand_file,    
                                lig_atom_ids=None, lig_atom_id_map=None):    
    """    
    为配体内部空间边分配相互作用类型    
    """    
    if len(lig_sei) == 0:    
        return lig_sei, lig_sea    
        
    try:     
            
        # 使用OpenBabel加载配体结构    
        molecule = pybel.readfile("sdf", ligand_file).__next__()    
            
        # 创建并设置分子，使用精细检测函数    
        analyzer = LigandInteractionAnalyzer()    
        analyzer.set_molecule(molecule)  # 预计算所有功能基团    
          
        # 方案1：预构建坐标到原子的哈希映射（高效且准确）  
        coord_to_atom_map = {}  
          
        for atom in molecule.atoms:  
            coord_key = tuple(np.round(atom.coords, 2))  # 保留2位小数精度  
            coord_to_atom_map[coord_key] = atom  
          
        def find_atom_by_coord_fast(coord):  
            coord_key = tuple(np.round(coord, 2))  
            return coord_to_atom_map.get(coord_key)  
            
        # 为每条空间边分配类型    
        classified_edge_attr = []  
        successful_matches = 0  
        failed_matches = 0  
        interaction_detected = 0  
        total_edges = len(lig_sei.T)  
            
        for i, (src_idx, tgt_idx) in enumerate(lig_sei.T):    
            src_coord = lig_coord[src_idx]    
            tgt_coord = lig_coord[tgt_idx]    
            distance = np.linalg.norm(src_coord - tgt_coord)    
                
            # 使用快速坐标查找  
            src_ob_atom = find_atom_by_coord_fast(src_coord)  
            tgt_ob_atom = find_atom_by_coord_fast(tgt_coord)  
              
            src_atom = src_ob_atom if src_ob_atom else None  
            tgt_atom = tgt_ob_atom if tgt_ob_atom else None  
                
            if src_atom and tgt_atom:  
                successful_matches += 1  
                # 使用精细的LigandInteractionAnalyzer分析原子对    
                interaction_result = analyzer.analyze_atom_pair(src_atom, tgt_atom, distance)    
                    
                if interaction_result and interaction_result['interaction_types']:  
                    interaction_detected += 1  
                    # 取第一个检测到的相互作用类型    
                    interaction_type = interaction_result['interaction_types'][0]    
                      
                    # if interaction_detected <= 10:  
                        # print(f"[DEBUG] 配体边 {i}: 检测到相互作用 {interaction_type}, 距离: {distance:.2f}Å")  
                        
                    # 映射到标准类型名称    
                    if interaction_type == 'hydrogen_bond':    
                        edge_type_name = 'hydrogen_bonds'    
                    elif interaction_type == 'hydrophobic':    
                        edge_type_name = 'hydrophobic_contacts'    
                    elif interaction_type == 'salt_bridge':    
                        edge_type_name = 'salt_bridges'    
                    elif interaction_type == 'halogen_bond':    
                        edge_type_name = 'halogen_bonds'    
                    else:    
                        edge_type_name = 'others'    
                else:    
                    # 检查π-π堆积和π-阳离子相互作用，使用预计算的功能基团    
                    edge_type_name = check_ring_interactions_improved(    
                        src_atom, tgt_atom, distance,     
                        analyzer.rings, analyzer.charged_groups    
                    )    
                        
                    # 检查水桥相互作用    
                    if edge_type_name == 'others' and detect_water_bridges_from_atoms(src_atom, tgt_atom, distance, molecule):    
                        edge_type_name = 'water_bridges'    
                        
                    # 检查金属配位相互作用    
                    if edge_type_name == 'others' and detect_metal_complex_from_atoms(src_atom, tgt_atom, distance, molecule):    
                        edge_type_name = 'metal_complexes'    
            else:  
                failed_matches += 1  
                edge_type_name = 'others'    
                
            # 获取边类型编码 - 使用当前文件中已定义的常量    
            INTERACTION_TYPE_MAP = {    
                'hydrogen_bonds': HYDROGEN_BOND_EDGE,    
                'hydrophobic_contacts': HYDROPHOBIC_EDGE,    
                'pi_stacking': PI_STACKING_EDGE,    
                'pi_cation': PI_CATION_EDGE,    
                'salt_bridges': SALT_BRIDGE_EDGE,    
                'water_bridges': WATER_BRIDGE_EDGE,    
                'halogen_bonds': HALOGEN_BOND_EDGE,    
                'metal_complexes': METAL_COMPLEX_EDGE,    
                'others': OTHERS_EDGE    
            }    
                
            edge_type = INTERACTION_TYPE_MAP.get(edge_type_name, OTHERS_EDGE)    
                
            # 创建边特征：[类型编码] + [距离]    
            edge_feature = edge_type + [distance]    
            classified_edge_attr.append(edge_feature)  
          
        # 添加统计信息  
        # print(f"[DEBUG] 配体原子匹配统计:")  
        # print(f"  成功匹配: {successful_matches}/{total_edges} ({successful_matches/total_edges*100:.1f}%)")  
        # print(f"  匹配失败: {failed_matches}/{total_edges} ({failed_matches/total_edges*100:.1f}%)")  
        # print(f"  检测到特定相互作用: {interaction_detected}/{successful_matches} ({interaction_detected/successful_matches*100:.1f}% of successful matches)")  
            
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)    
            
    except Exception as e:    
        print(f"配体相互作用分析失败: {e}")    
        # 如果分析失败，将所有边设置为OTHERS类型    
        classified_edge_attr = []    
        for i, (src_idx, tgt_idx) in enumerate(lig_sei.T):    
            distance = np.linalg.norm(lig_coord[src_idx] - lig_coord[tgt_idx])    
            edge_feature = OTHERS_EDGE + [distance]    
            classified_edge_attr.append(edge_feature)    
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)    
        
    return lig_sei, classified_edge_attr

def find_closest_atom_by_coord(coord, molecule):  
    """根据坐标找到最接近的原子"""  
    min_dist = float('inf')  
    closest_atom = None  
      
    for atom in molecule.atoms:  
        atom_coord = np.array(atom.coords)  
        dist = np.linalg.norm(np.array(coord) - atom_coord)  
        if dist < min_dist:  
            min_dist = dist  
            closest_atom = AtomInfo(atom)  
      
    return closest_atom if min_dist < 2.0 else None  # 2Å阈值  
  
def check_ring_interactions_improved(src_atom, tgt_atom, distance, rings, charged_groups):  
    """改进的环相互作用检测，使用PLIP的精细检测函数"""  
    from intra_pro_plip import pistacking, pication  
      
    # 找到包含这两个原子的环  
    src_rings = []  
    tgt_rings = []  
      
    for ring in rings:  
        if any(ring_atom.idx == src_atom.idx for ring_atom in ring.atoms):  
            src_rings.append(ring)  
        if any(ring_atom.idx == tgt_atom.idx for ring_atom in ring.atoms):  
            tgt_rings.append(ring)  
      
    # 使用PLIP的精细π-π堆积检测  
    if src_rings and tgt_rings:  
        pi_stack_results = pistacking(src_rings, tgt_rings)  
        if pi_stack_results:  
            return 'pi_stacking'  
      
    # 使用PLIP的精细π-阳离子检测  
    if distance < 6.0:  # PICATION_DIST_MAX  
        # 检查是否有环-电荷组合  
        for ring in src_rings + tgt_rings:  
            for charge_group in charged_groups:  
                if charge_group.type == 'positive':  
                    pi_cation_results = pication([ring], [charge_group])  
                    if pi_cation_results:  
                        return 'pi_cation'  
      
    return 'others'

def set_ligand_spatial_edge_types(lig_sea):  
    """  
    将配体内部空间边的类型统一设置为OTHERS_EDGE  
      
    Args:  
        lig_sea: 配体空间边属性  
      
    Returns:  
        np.array: 更新后的边属性  
    """  
    if len(lig_sea) == 0:  
        return lig_sea  
      
    # 为每条边设置OTHERS类型  
    updated_edge_attr = []  
      
    for edge_attr in lig_sea:  
        # 保留原有的距离信息（如果存在）  
        if len(edge_attr) > 3:  
            distance = edge_attr[3]  
        else:  
            distance = edge_attr[-1] if len(edge_attr) > 0 else 0.0  
          
        # 创建新的边特征：OTHERS_EDGE + 距离  
        new_edge_feature = OTHERS_EDGE + [distance]  
        updated_edge_attr.append(new_edge_feature)  
      
    return np.array(updated_edge_attr, dtype=np.float32)

def gen_graph(ligand: tuple, pocket: tuple, name: str, protein_cutoff: float,     
              pocket_cutoff: float, spatial_cutoff: float, protein_file=None,     
              ligand_file=None, plip_interactions=None,  
              ligand_dict=None, pocket_dict=None):  # 新增参数   
    # 构建复合物的整体图结构（节点、特征、边、空间边）  
  
    # 解包配体和蛋白 pocket 的输入：坐标、特征、边索引、边特征  
    lig_coord, lig_feat, lig_ei, lig_ea = ligand  
    pro_coord, pro_feat, pro_ei, pro_ea = pocket  
  
    # 检查配体和蛋白的坐标和特征长度是否一致  
    assert len(lig_coord) == len(lig_feat)  
    assert len(pro_coord) == len(pro_feat)  
  
    # 检查 cutoff 参数的合理性  
    assert protein_cutoff >= pocket_cutoff, f"Protein cutoff {protein_cutoff} should be larger than pocket cutoff {pocket_cutoff}"  
    assert pocket_cutoff >= spatial_cutoff, f"Pocket cutoff {pocket_cutoff} should be larger than spatial cutoff {spatial_cutoff}"  
  
    # 如果没有提供PLIP分析结果且有文件路径，进行PLIP分析  
    if plip_interactions is None and protein_file and ligand_file:  
        plip_interactions = analyze_plip_interactions(protein_file, ligand_file)  
    # 只保留距离配体小于 protein_cutoff 的蛋白原子
    # 先初始化全为 False 的掩码
    pro_atom_mask = np.zeros(len(pro_coord), dtype=bool)
    # 计算配体和蛋白的距离矩阵，保留距离小于 protein_cutoff 的蛋白原子的索引
    pro_atom_mask[np.where(cdist(lig_coord, pro_coord) <= protein_cutoff)[1]] = 1
    # 只保留蛋白内部边的两端都在掩码范围内的边
    pro_edge_mask = np.array([True if pro_atom_mask[i] and pro_atom_mask[j] else False for i, j in zip(*pro_ei)])

    # 筛选对应的蛋白节点和边
    pro_coord = pro_coord[pro_atom_mask]   # 只保留掩码内的坐标
    pro_feat = pro_feat[pro_atom_mask]     # 只保留掩码内的特征
    # 只保留掩码内的边
    pro_ei = np.array([pro_ei[0, pro_edge_mask], pro_ei[1, pro_edge_mask]])
    pro_ea = pro_ea[pro_edge_mask]
    # 边的原子索引需要重新编号（因为节点删减了）
    pro_ei = reindex(np.where(pro_atom_mask)[0], pro_ei)
    lig_sei = np.empty((2, 0), dtype=np.int64)  
    lig_sea = np.empty((0, 4), dtype=np.float32)
    pro_sei = np.empty((2, 0), dtype=np.int64)  
    pro_sea = np.empty((0, 4), dtype=np.float32)
    # 配体内部空间边：所有距离小于 spatial_cutoff 的配体原子对
    lig_dm = cdist(lig_coord, lig_coord)   # 配体原子之间的距离矩阵
    lig_atom_ids = ligand_dict.get('atom_ids') if ligand_dict else None  
    lig_atom_id_map = ligand_dict.get('atom_id_to_graph_index') if ligand_dict else None  
    pro_atom_ids = pocket_dict.get('atom_ids') if pocket_dict else None  
    pro_atom_id_map = pocket_dict.get('atom_id_to_graph_index') if pocket_dict else None 
    lig_sei, lig_sea = gen_spatial_edge(lig_dm, spatial_cutoff=spatial_cutoff)  # 生成空间边
    lig_sei, lig_sea = remove_duplicated_edges(lig_sei, lig_sea, lig_ei)        # 去掉和结构边重复的空间边
    lig_sei, lig_sea = classify_ligand_spatial_edges(  
        lig_sei, lig_sea, lig_coord, ligand_file,  
        lig_atom_ids=lig_atom_ids, lig_atom_id_map=lig_atom_id_map  
    ) 
    # 蛋白内部空间边：所有距离小于 spatial_cutoff 的蛋白原子对
    pro_dm = cdist(pro_coord, pro_coord)  
    pro_sei, pro_sea = gen_spatial_edge(pro_dm, spatial_cutoff=spatial_cutoff)  
    pro_sei, pro_sea = remove_duplicated_edges(pro_sei, pro_sea, pro_ei)
    pro_sei, pro_sea = classify_protein_spatial_edges(  
        pro_sei, pro_sea, pro_coord, protein_file,  
        pro_atom_ids=pro_atom_ids, pro_atom_id_map=pro_atom_id_map  
    )  
    # 新增：使用更新的 intra_pro_plip 脚本进行类型分类  
    
    # 配体-蛋白之间的空间边：所有距离小于 pocket_cutoff 的配体-蛋白原子对
    dm_lig_pro = cdist(lig_coord, pro_coord)  
    lig_pock_ei, lig_pock_ea = gen_ligpro_edge(  
        dm_lig_pro,  
        pocket_cutoff=pocket_cutoff,  
        plip_interactions=plip_interactions,  
        lig_coord=lig_coord,  
        pro_coord=pro_coord,  
        lig_atom_id_map=lig_atom_id_map,  # 新增  
        pro_atom_id_map=pro_atom_id_map   # 新增  
    ) 
    def pad_edge_features_with_distance(edge_attr, edge_index, coords):  
        """为3维边特征添加距离信息，使其变为4维"""  
        if edge_attr.shape[1] == 3:  # 如果是3维特征  
            distances = []  
            for i in range(edge_index.shape[1]):  
                src_idx, dst_idx = edge_index[0, i], edge_index[1, i]  
                dist = np.linalg.norm(coords[src_idx] - coords[dst_idx])  
                distances.append(dist)  

            # 添加距离作为第4维  
            distances = np.array(distances).reshape(-1, 1)  
            edge_attr = np.hstack([edge_attr, distances]).astype(np.float32)  

        return edge_attr  
  
# 在拼接前处理所有边特征  
    lig_ea = pad_edge_features_with_distance(lig_ea, lig_ei, lig_coord)  
    pro_ea = pad_edge_features_with_distance(pro_ea, pro_ei, pro_coord)
    # 构建整体节点、特征
    comp_coord = np.vstack([lig_coord, pro_coord])   # 合并配体和蛋白的原子坐标
    comp_feat = np.vstack([lig_feat, pro_feat])      # 合并配体和蛋白的原子特征
    comp_ei, comp_ea = lig_ei, lig_ea                # 先只用配体的结构边作为初始边
    # print(f"lig_ea shape: {lig_ea.shape}")  
    # print(f"pro_ea shape: {pro_ea.shape}")  
    # print(f"lig_sea shape: {lig_sea.shape}")  
    # print(f"pro_sea shape: {pro_sea.shape}")  
    # print(f"lig_pock_ea shape: {lig_pock_ea.shape}")

    # 拼接所有蛋白、配体、空间和配体-蛋白边
    if len(pro_ei.shape) == 2 and len(pro_ei.T) >= 3:   # 如果蛋白有结构边
        comp_ei = np.hstack([comp_ei, pro_ei + len(lig_feat)])   # 蛋白边下标整体平移
        comp_ea = np.vstack([comp_ea, pro_ea])                  # 合并蛋白边特征
    if len(lig_sei.shape) == 2 and len(lig_sei.T) >= 3:         # 如果配体有空间边
        comp_ei = np.hstack([comp_ei, lig_sei])                 # 合并配体空间边
        comp_ea = np.vstack([comp_ea, lig_sea])                 # 合并配体空间边特征
    if len(pro_sei.shape) == 2 and len(pro_sei.T) >= 3:         # 如果蛋白有空间边
        comp_ei = np.hstack([comp_ei, pro_sei + len(lig_feat)]) # 蛋白空间边下标整体平移
        comp_ea = np.vstack([comp_ea, pro_sea])                 # 合并蛋白空间边特征
    if len(lig_pock_ei.shape) == 2 and len(lig_pock_ei.T) >= 3: # 如果有配体-蛋白空间边
        comp_ei = np.hstack([comp_ei, lig_pock_ei])             # 合并配体-蛋白空间边
        comp_ea = np.vstack([comp_ea, lig_pock_ea])             # 合并配体-蛋白空间边特征

    # 记录每类节点、边数
    comp_num_node = np.array([len(lig_feat), len(pro_feat)], dtype=np.int64)   # 配体、蛋白节点数
    comp_num_edge = np.array([  
    lig_ei.shape[1],           # 配体结构边数  
    pro_ei.shape[1],           # 蛋白结构边数  
    lig_pock_ei.shape[1],      # 配体-蛋白边数  
    lig_sei.shape[1],          # 配体空间边数
    pro_sei.shape[1],          # 蛋白空间边数
], dtype=np.int64)
    return comp_coord, comp_feat, comp_ei, comp_ea, comp_num_node, comp_num_edge, lig_sei, lig_sea, pro_sei, pro_sea

def load_pk_data(data_path: Path):
    # 从txt读取pdbid与pk值
    pdbid, pk = [], []
    for line in open(data_path):
        if line[0] == '#': continue
        elem = line.split()
        v1, _, _, v2 = elem[:4]
        pdbid.append(v1)
        pk.append(float(v2))
    res = {i: p for i, p in zip(pdbid, pk)}
    return res

# def to_pyg_graph(raw: list, **kwargs):  
#     comp_coord, comp_feat, comp_ei, comp_ea, comp_num_node, comp_num_edge, rfscore, gbscore, ecif, pk, name = raw  
      
#     # 从kwargs中提取rmsd值  
#     rmsd = kwargs.pop('rmsd', 0.0)  # 默认值为0.0  
  
#     d = Data(x=torch.from_numpy(comp_feat).to(torch.long),     
#          edge_index=torch.from_numpy(comp_ei).to(torch.long),     
#          edge_attr=torch.from_numpy(comp_ea).to(torch.float32),  # 改为float32  
#          pos=torch.from_numpy(comp_coord).to(torch.float32),
#      y=torch.tensor([pk], dtype=torch.float32),   
#      pdbid=name,  
#      num_node=torch.from_numpy(comp_num_node).to(torch.long),   
#      num_edge=torch.from_numpy(comp_num_edge).to(torch.long),  
#      rfscore=torch.from_numpy(rfscore).to(torch.float32),   
#      gbscore=torch.from_numpy(gbscore).to(torch.float32),  
#      ecif=torch.from_numpy(ecif).to(dtype=torch.float32),  
#      rmsd=torch.tensor([rmsd], dtype=torch.float32),  # 新增：rmsd值  
#      **kwargs)  
#     return d

def get_info(protein_file, ligand_file):
    # 用PyMOL获取蛋白和配体的元素、残基和坐标信息
    cmd.reinitialize()
    cmd.load(protein_file, 'receptor')
    cmd.load(ligand_file, 'ligand')
    cmd.remove('sol.')
    cmd.h_add()
    proinfo = {"elem": [], "resn": [], "coord":[], }
    liginfo = {"elem": [], "resn": [], "coord":[], }
    cmd.iterate_state(1, 'receptor', 'info["elem"].append(elem); info["resn"].append(resn); info["coord"].append(np.array([x, y, z]))',space={"info": proinfo, "np": np})
    cmd.iterate_state(1, 'ligand', 'info["elem"].append(elem); info["resn"].append(resn); info["coord"].append(np.array([x, y, z]))',space={"info": liginfo, "np": np})
    for k in proinfo.keys():
        proinfo[k] = np.array(proinfo[k])
        liginfo[k] = np.array(liginfo[k])
    return proinfo, liginfo

def GB_score(lig_info: dict, pro_info: dict) -> np.ndarray:
    # 生成GBScore特征（配体-蛋白原子对分组距离统计，400维）
    amino_acid_groups = [
        {"ARG", "LYS", "ASP", "GLU"},
        {"GLN", "ASN", "HIS", "SER", "THR", "CYS"},
        {"TRP", "TYR", "MET"},
        {"ILE", "LEU", "PHE", "VAL", "PRO", "GLY", "ALA"},
    ]
    elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
    distmap = cdist(lig_info['coord'], pro_info['coord'])
    restype = np.zeros(len(pro_info['resn'])) - 1
    elem_mask = {k: pro_info["elem"] == k for k in elements}
    fp = np.zeros([len(elements), len(elements), len(amino_acid_groups)])

    for idx, r in enumerate(pro_info['resn']):
        if r in amino_acid_groups[0]: restype[idx] = 0
        if r in amino_acid_groups[1]: restype[idx] = 1
        if r in amino_acid_groups[2]: restype[idx] = 2
        if r in amino_acid_groups[3]: restype[idx] = 3

    for i, el in enumerate(elements):
        lmask = lig_info["elem"] == el
        if lmask.sum() < 1: continue
        for j, ep in enumerate(elements):
            pmask = elem_mask[ep]
            if pmask.sum() < 1: continue
            for k, rt in enumerate(range(4)):
                rt_mask = restype == rt
                m = pmask & rt_mask
                if m.sum() < 1: continue
                d = distmap[lmask][:, m]
                v = (1 / d[d<=12]).sum()
                fp[i, j, k] = v

    return fp.flatten()

