# 从 features.py 文件中导入若干函数和变量，用于分子特征处理
from features import (
    allowable_features,             # 允许的特征列表
    atom_to_feature_vector,         # 原子对象转特征向量
    bond_to_feature_vector,         # 键对象转特征向量
    atom_feature_vector_to_dict,    # 原子特征向量转字典
    bond_feature_vector_to_dict     # 键特征向量转字典
)
from pathlib import Path
from rdkit import Chem
from utils import correct_sanitize_v2                # 导入 RDKit 的 Chem 模块（化学结构处理）
import numpy as np                 # 导入 numpy 用于数值计算
from utils import read_mol         # 从 utils.py 中导入 read_mol 函数
from features import atom_to_feature_vector, bond_to_feature_vector  

def mol2graph(mol: Chem.Mol):
    conformer = mol.GetConformer(0)    # 获取分子的第一个构象对象，用于获取3D坐标

    # 处理原子
    atom_features_list, coords = [], []    # 用于存储原子特征和原子坐标
    atom_map = dict()                      # 建立原子在新序号中的映射（去除氢原子）
    for idx, atom in enumerate(mol.GetAtoms()):   # 遍历分子中的每个原子
        if atom.GetSymbol() == "H":        # 如果是氢原子则跳过
            continue
        atom_features_list.append(atom_to_feature_vector(atom))  # 原子特征向量化
        coords.append(conformer.GetAtomPosition(atom.GetIdx()))  # 获取原子坐标
        atom_map[idx] = len(coords) - 1    # 原子在新表中的索引（跳过氢原子）

    x = np.array(atom_features_list, dtype = np.int64)   # 所有原子特征转换为 numpy 数组

    # 处理化学键
    num_bond_features = 3  # 键的特征数量（类型、立体、是否共轭）
    if len(mol.GetBonds()) > 0:    # 如果分子有化学键
        edges_list = []            # 存储边（键）两端的原子索引
        edge_features_list = []    # 存储边的特征
        for bond in mol.GetBonds():    # 遍历分子中的每个化学键
            i = bond.GetBeginAtomIdx()   # 键的起始原子索引
            j = bond.GetEndAtomIdx()     # 键的终止原子索引
            # 如果任一端是氢原子则跳过
            if mol.GetAtomWithIdx(i).GetSymbol() == "H":
                continue
            if mol.GetAtomWithIdx(j).GetSymbol() == "H":
                continue

            edge_feature = bond_to_feature_vector(bond)    # 键特征向量化

            # 添加有向边：i->j 和 j->i（无向图）
            edges_list.append((atom_map[i], atom_map[j]))
            edge_features_list.append(edge_feature)
            edges_list.append((atom_map[j], atom_map[i]))
            edge_features_list.append(edge_feature)

        # 生成 numpy 数组，shape=[2,num_edges]，表示 COO 格式的边索引
        edge_index = np.array(edges_list, dtype = np.int64).T
        # 边的特征，shape=[num_edges, num_edge_features]
        edge_attr = np.array(edge_features_list, dtype = np.int64)

    else:   # 如果分子没有化学键
        edge_index = np.empty((2, 0), dtype = np.int64)        # 空边索引
        edge_attr = np.empty((0, num_bond_features), dtype = np.int64)  # 空边特征

    # 构建分子图字典
    graph = dict()
    graph['edge_index'] = edge_index       # 边索引
    graph['edge_feat'] = edge_attr         # 边特征
    graph['node_feat'] = x                 # 节点（原子）特征
    graph['coords'] = np.array(coords)     # 原子坐标

    return graph  # 返回分子图

def mol2graph_protein_from_pdb(pdb_file: Path):  
    """直接从PDB文件处理蛋白质分子的图构建"""  
    # 使用RDKit读取PDB文件  
    mol = Chem.MolFromPDBFile(str(pdb_file), sanitize=False, removeHs=False)  
    if mol is None:  
        raise RuntimeError(f"Cannot read PDB file: {pdb_file}")  
      
    mol = correct_sanitize_v2(mol)  
    mol = Chem.RemoveHs(mol, sanitize=False)  
      
    conformer = mol.GetConformer(0)  
    atom_features_list, coords, pro_names, aa_names = [], [], [], []  
    atom_map = dict()  
      
    for idx, atom in enumerate(mol.GetAtoms()):  
        if atom.GetSymbol() == "H": continue  
          
        atom_features_list.append(atom_to_feature_vector(atom))  
        coords.append(conformer.GetAtomPosition(atom.GetIdx()))  
          
        # 从PDB残基信息中获取原子名称和氨基酸名称  
        atom_info = atom.GetPDBResidueInfo()  
        if atom_info:  
            pro_names.append(atom_info.GetName().strip())  
            aa_names.append(atom_info.GetResidueName().strip())  
        else:  
            pro_names.append("UNK")  
            aa_names.append("UNK")  
              
        atom_map[idx] = len(coords) - 1  
      
    x = np.array(atom_features_list, dtype=np.int64)  
      
    # 处理边信息（完整的边处理逻辑）  
    num_bond_features = 3  # bond type, bond stereo, is_conjugated  
    if len(mol.GetBonds()) > 0: # mol has bonds  
        edges_list = []  
        edge_features_list = []  
        for bond in mol.GetBonds():  
            i = bond.GetBeginAtomIdx()  
            j = bond.GetEndAtomIdx()  
            if mol.GetAtomWithIdx(i).GetSymbol() == "H": continue  
            if mol.GetAtomWithIdx(j).GetSymbol() == "H": continue  
  
            edge_feature = bond_to_feature_vector(bond)  
  
            # add edges in both directions  
            edges_list.append((atom_map[i], atom_map[j]))  
            edge_features_list.append(edge_feature)  
            edges_list.append((atom_map[j], atom_map[i]))  
            edge_features_list.append(edge_feature)  
  
        # data.edge_index: Graph connectivity in COO format with shape [2, num_edges]  
        edge_index = np.array(edges_list, dtype=np.int64).T  
  
        # data.edge_attr: Edge feature matrix with shape [num_edges, num_edge_features]  
        edge_attr = np.array(edge_features_list, dtype=np.int64)  
  
    else:   # mol has no bonds  
        edge_index = np.empty((2, 0), dtype=np.int64)  
        edge_attr = np.empty((0, num_bond_features), dtype=np.int64)  
      
    graph = dict()  
    graph['edge_index'] = edge_index  
    graph['edge_feat'] = edge_attr  
    graph['node_feat'] = x  
    graph['coords'] = np.array(coords)  
    graph['pro_name'] = np.array(pro_names)  
    graph['AA_name'] = np.array(aa_names)  
      
    return graph
  
def mol2graph_ligand(mol: Chem.Mol):  
    """处理小分子的图构建"""  
    # 获取SMILES  
    smiles = Chem.MolToSmiles(mol)  
      
    # 调用原有逻辑  
    result = mol2graph(mol)  
    result['smiles'] = smiles  # 新增：SMILES字符串  
      
    return result  
# 主程序入口
if __name__ == '__main__':
    graph = smiles2graph('O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5')
    print(graph)    # 输出分子图