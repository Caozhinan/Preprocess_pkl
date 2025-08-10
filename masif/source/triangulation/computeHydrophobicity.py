import numpy as np  
from rdkit import Chem  
from rdkit.Chem import Crippen  
from rdkit.Chem import BRICS  
  
  
# Kyte Doolittle scale  
kd_scale = {}  
kd_scale["ILE"] = 4.5  
kd_scale["VAL"] = 4.2  
kd_scale["LEU"] = 3.8  
kd_scale["PHE"] = 2.8  
kd_scale["CYS"] = 2.5  
kd_scale["MET"] = 1.9  
kd_scale["ALA"] = 1.8  
kd_scale["GLY"] = -0.4  
kd_scale["THR"] = -0.7  
kd_scale["SER"] = -0.8  
kd_scale["TRP"] = -0.9  
kd_scale["TYR"] = -1.3  
kd_scale["PRO"] = -1.6  
kd_scale["HIS"] = -3.2  
kd_scale["GLU"] = -3.5  
kd_scale["GLN"] = -3.5  
kd_scale["ASP"] = -3.5  
kd_scale["ASN"] = -3.5  
kd_scale["LYS"] = -3.9  
kd_scale["ARG"] = -4.5  
  
  
def kd_from_logp(logp, kd_min=-4.5, kd_max=4.5):  
    return np.clip(-6.2786 + np.exp(0.4772 * logp + 1.8491), kd_min, kd_max)  
  
  
def get_fragments(rdmol):  
    frags = Chem.GetMolFrags(Chem.FragmentOnBRICSBonds(rdmol), asMols=True, sanitizeFrags=False)  
  
    # find exit atoms to be removed  
    def find_exits(rdmol):  
        exits = []  
        for a in rdmol.GetAtoms():  
            if a.GetSymbol() == '*':  
                exits.append(a.GetIdx())  
        return exits  
  
    out_frags = []  
    for frag in frags:  
        exits = sorted(find_exits(frag), reverse=True)  
        efrag = Chem.EditableMol(frag)  
        for idx in exits:  
            efrag.RemoveAtom(idx)  
        efrag = efrag.GetMol()  
        Chem.SanitizeMol(efrag, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)  
        out_frags.append(efrag)  
    return out_frags  
  
  
def custom_mol_assign_kd(mol, msms_atom_names):  
    """  
    自定义的分子疏水性分配函数  
    直接基于MSMS原子名称进行映射，避免PDB残基信息依赖  
    """  
    # 计算整个分子的logP  
    mol_logp = Chem.Crippen.MolLogP(mol)  
    mol_kd = kd_from_logp(mol_logp)  
      
    # 为所有MSMS原子名称分配相同的疏水性值  
    atom_kd = {}  
    for atom_name in set(msms_atom_names):  
        atom_kd[atom_name] = mol_kd  
      
    return atom_kd  
  
  
# For each vertex in names, compute  
def computeHydrophobicity(names, ligand_code=None, rdmol=None):  
    """  
    修改后的疏水性计算函数，避免原子名称映射问题  
    """  
    hp = np.zeros(len(names))  
      
    if rdmol is not None:  
        # 提取所有配体原子的MSMS名称  
        ligand_atom_names = []  
        for ix, name in enumerate(names):  
            aa = name.split("_")[3]  
            if aa == ligand_code:  
                atom_name = name.split("_")[4]  
                ligand_atom_names.append(atom_name)  
          
        # 使用自定义映射函数  
        mol_kd = custom_mol_assign_kd(rdmol, ligand_atom_names)  
      
    for ix, name in enumerate(names):  
        aa = name.split("_")[3]  
          
        if rdmol is not None and aa == ligand_code:  
            atom_name = name.split("_")[4]  
            hp[ix] = mol_kd.get(atom_name, 0.0)  # 使用get避免KeyError  
        else:  
            # 为未知残基类型添加默认值  
            hp[ix] = kd_scale.get(aa, 0.0)
      
    return hp