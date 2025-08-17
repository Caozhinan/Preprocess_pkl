#!/usr/bin/env python3    
"""    
配体内部相互作用分析模块    
直接复用 intra_pro_plip.py 的所有类和函数定义，适配配体分子    
"""    
    
import numpy as np    
from openbabel import pybel    
from typing import List, Dict, Tuple, Optional, Set    
import math  
import itertools  
from collections import namedtuple  
  
# 导入 intra_pro_plip.py 中的精细检测函数和辅助函数  
from intra_pro_plip import (  
    hydrophobic_interactions, hbonds, pistacking,   
    saltbridge, pication, halogen_bonds,  
    euclidean3d, vector, vecangle, normalize_vector,  
    AtomInfo as ProAtomInfo, RingInfo as ProRingInfo  
)  
  
class AtomInfo:    
    """原子信息类，存储原子的基本属性 - 兼容 intra_pro_plip.py"""    
    def __init__(self, atom):    
        self.atom = atom    
        self.idx = atom.idx    
        self.atomicnum = atom.atomicnum    
        self.coords = np.array([atom.coords[0], atom.coords[1], atom.coords[2]])    
        self.element = atom.OBAtom.GetAtomicNum()    
        self.formal_charge = getattr(atom, 'formalcharge', 0)
        self.partial_charge = getattr(atom.OBAtom, 'GetPartialCharge', lambda: 0.0)()  
        self.type = atom.type  
        self.OBAtom = atom.OBAtom  
  
class Config:    
    """配置参数类 - 与 intra_pro_plip.py 保持一致"""    
    HYDROPHOBIC_DIST_MAX = 4.0  
    HYDROPH_DIST_MAX = 4.0  # 别名，与 intra_pro_plip.py 兼容  
    HBOND_DIST_MAX = 3.5    
    HBOND_DON_ANGLE_MIN = 120.0    
    PISTACK_DIST_MAX = 5.5    
    PISTACK_ANG_DEV = 30.0    
    PISTACK_OFFSET_MAX = 2.0    
    PICATION_DIST_MAX = 6.0    
    SALTBRIDGE_DIST_MAX = 5.0    
    HALOGEN_DIST_MAX = 4.0    
    HALOGEN_ACC_ANGLE = 120.0    
    HALOGEN_DON_ANGLE = 165.0  
    HALOGEN_ANGLE_DEV = 30.0  # 新增角度偏差参数  
    WATER_BRIDGE_MINDIST = 2.5    
    WATER_BRIDGE_MAXDIST = 4.0    
    METAL_DIST_MAX = 3.0  
    MIN_DIST = 0.5  # 最小距离阈值  
  
# 全局配置实例  
config = Config()  
  
def find_hydrophobic_atoms(molecule):    
    """识别疏水原子 - 返回 AtomInfo 对象列表，兼容 intra_pro_plip.py"""    
    hydrophobic_atoms = []    
        
    for atom in molecule.atoms:    
        if atom.atomicnum == 6:  # Carbon    
            # 检查邻居原子是否只有碳和氢    
            neighbor_nums = set()    
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):    
                neighbor_nums.add(neighbor.GetAtomicNum())    
                
            # 如果邻居只有碳(6)和氢(1)，则为疏水原子    
            if neighbor_nums.issubset({1, 6}):    
                hydrophobic_atoms.append(AtomInfo(atom))    
        elif atom.atomicnum == 16:  # Sulfur    
            hydrophobic_atoms.append(AtomInfo(atom))    
        
    return hydrophobic_atoms    
  
def find_hba(molecule):    
    """识别氢键受体 - 返回 AtomInfo 对象列表，兼容 intra_pro_plip.py"""    
    acceptors = []    
        
    for atom in molecule.atoms:    
        # 排除卤素原子    
        if atom.atomicnum not in [9, 17, 35, 53]:    
            if atom.OBAtom.IsHbondAcceptor():    
                acceptors.append(AtomInfo(atom))    
        
    return acceptors    
  
def find_hbd(molecule, hydrophobic_atoms=None):    
    """识别氢键供体 - 返回 namedtuple 对象列表，兼容 intra_pro_plip.py"""    
    data = namedtuple('hbonddonor', 'd h type')    
    donors = []    
        
    # 强氢键供体    
    for atom in molecule.atoms:    
        if atom.OBAtom.IsHbondDonor():    
            for adj_atom in pybel.ob.OBAtomAtomIter(atom.OBAtom):    
                if adj_atom.IsHbondDonorH():    
                    donors.append(data(d=AtomInfo(atom), h=pybel.Atom(adj_atom), type='regular'))    
        
    # 弱氢键供体（疏水碳-氢）    
    if hydrophobic_atoms:    
        for carbon in hydrophobic_atoms:    
            for adj_atom in pybel.ob.OBAtomAtomIter(carbon.atom.OBAtom):    
                if adj_atom.GetAtomicNum() == 1:  # 氢原子    
                    donors.append(data(d=carbon, h=pybel.Atom(adj_atom), type='weak'))    
        
    return donors    
  
class RingInfo:    
    """芳香环信息类 - 兼容 intra_pro_plip.py"""    
    def __init__(self, atoms, center, normal, ring_type='aromatic'):    
        self.atoms = atoms  # AtomInfo 对象列表    
        self.center = center    
        self.normal = normal    
        self.type = ring_type    
  
def find_rings(molecule):    
    """识别芳香环 - 返回 RingInfo 对象列表，兼容 intra_pro_plip.py"""    
    rings = []    
        
    ring_candidates = molecule.OBMol.GetSSSR()    
        
    for ring in ring_candidates:    
        ring_atoms = []    
        for atom in molecule.atoms:    
            if ring.IsMember(atom.OBAtom):    
                ring_atoms.append(AtomInfo(atom))    
            
        ring_atoms = sorted(ring_atoms, key=lambda x: x.idx)    
            
        if 4 < len(ring_atoms) <= 6:    
            if ring.IsAromatic():    
                ring_coords = [atom.coords for atom in ring_atoms]    
                center = np.mean(ring_coords, axis=0)    
                    
                if len(ring_atoms) >= 3:    
                    v1 = vector(ring_coords[0], ring_coords[2])    
                    v2 = vector(ring_coords[4 if len(ring_atoms) > 4 else 1], ring_coords[0])    
                    normal = normalize_vector(np.cross(v1, v2))    
                else:    
                    normal = np.array([0, 0, 1])    
                    
                ring_obj = RingInfo(ring_atoms, center, normal)    
                rings.append(ring_obj)    
        
    return rings    
  
def find_charged_groups(molecule):    
    """识别带电基团 - 返回 namedtuple 对象列表，兼容 intra_pro_plip.py"""    
    data = namedtuple('charge', 'atoms type center')    
    charged_groups = []    
        
    # 正电荷原子组    
    positive_atoms = []    
    negative_atoms = []    
        
    for atom in molecule.atoms:    
        if atom.formalcharge > 0:    
            positive_atoms.append(AtomInfo(atom))    
        elif atom.formalcharge < 0:    
            negative_atoms.append(AtomInfo(atom))    
            
        # 四价氮（季铵）    
        if atom.atomicnum == 7:    
            neighbors = list(pybel.ob.OBAtomAtomIter(atom.OBAtom))    
            if len(neighbors) == 4:    
                positive_atoms.append(AtomInfo(atom))    
            
        # 羧基氧    
        elif atom.atomicnum == 8:    
            neighbors = list(pybel.ob.OBAtomAtomIter(atom.OBAtom))    
            for neighbor in neighbors:    
                if neighbor.GetAtomicNum() == 6:    
                    carbon_neighbors = list(pybel.ob.OBAtomAtomIter(neighbor))    
                    oxygen_count = sum(1 for n in carbon_neighbors if n.GetAtomicNum() == 8)    
                    if oxygen_count >= 2:    
                        negative_atoms.append(AtomInfo(atom))    
                        break    
      
    # 构建带电基团对象    
    if positive_atoms:    
        center = np.mean([atom.coords for atom in positive_atoms], axis=0)    
        charged_groups.append(data(atoms=positive_atoms, type='positive', center=center))    
          
    if negative_atoms:    
        center = np.mean([atom.coords for atom in negative_atoms], axis=0)    
        charged_groups.append(data(atoms=negative_atoms, type='negative', center=center))    
        
    return charged_groups    
  
def find_halogen_acceptors(molecule):    
    """识别卤键受体 - 返回 namedtuple 对象列表，兼容 intra_pro_plip.py"""    
    data = namedtuple('hal_acceptor', 'o y')    
    acceptors = []    
        
    for atom in molecule.atoms:    
        if atom.atomicnum in [8, 7, 16]:  # O, N, S    
            neighbors = []    
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):    
                if neighbor.GetAtomicNum() in [6, 7, 15, 16]:  # C, N, P, S    
                    neighbors.append(pybel.Atom(neighbor))    
                
            if len(neighbors) >= 1:  # 至少有一个邻接原子    
                acceptors.append(data(o=AtomInfo(atom), y=neighbors[0]))    
        
    return acceptors    
  
def find_halogen_donors(molecule):    
    """识别卤键供体 - 返回 namedtuple 对象列表，兼容 intra_pro_plip.py"""    
    data = namedtuple('hal_donor', 'x c')    
    donors = []    
        
    for atom in molecule.atoms:    
        if atom.atomicnum in [9, 17, 35, 53]:  # F, Cl, Br, I    
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):    
                if neighbor.GetAtomicNum() == 6:  # 连接到碳原子    
                    donors.append(data(x=AtomInfo(atom), c=pybel.Atom(neighbor)))    
                    break    
        
    return donors    
  
def find_closest_atom_by_coord(coord: np.ndarray, molecule) -> Optional:     # type: ignore
    """根据坐标找到最接近的原子"""    
    min_dist = float('inf')    
    closest_atom = None    
        
    for atom in molecule.atoms:    
        atom_coord = np.array([atom.coords[0], atom.coords[1], atom.coords[2]])    
        dist = np.linalg.norm(coord - atom_coord)    
        if dist < min_dist:    
            min_dist = dist    
            closest_atom = atom    
        
    return closest_atom    
  
class InteractionAnalyzer:    
    """相互作用分析器 - 使用 intra_pro_plip.py 的精细检测函数"""    
        
    def __init__(self):    
        self.config = Config()  
        self.molecule = None  
        # 预计算的功能基团  
        self.hydrophobic_atoms = None  
        self.hba_atoms = None  
        self.hbd_atoms = None  
        self.rings = None  
        self.charged_groups = None  
        self.hal_acceptors = None  
        self.hal_donors = None  
      
    def set_molecule(self, molecule):  
        """设置分子并预计算所有功能基团"""  
        self.molecule = molecule  
        self.hydrophobic_atoms = find_hydrophobic_atoms(molecule)  
        self.hba_atoms = find_hba(molecule)  
        self.hbd_atoms = find_hbd(molecule, self.hydrophobic_atoms)  
        self.rings = find_rings(molecule)  
        self.charged_groups = find_charged_groups(molecule)  
        self.hal_acceptors = find_halogen_acceptors(molecule)  
        self.hal_donors = find_halogen_donors(molecule)  
        
    def analyze_atom_pair(self, atom1, atom2, distance: float) -> Optional[Dict]:    
        """分析原子对之间的相互作用 - 使用精细检测函数"""    
        if not self.molecule:  
            # 如果没有设置分子，回退到简化检测  
            return self._analyze_atom_pair_simple(atom1, atom2, distance)  
          
        interaction_types = []  
          
        # 转换为 AtomInfo 对象  
        atom1_info = AtomInfo(atom1)  
        atom2_info = AtomInfo(atom2)
        if self.hba_atoms and self.hbd_atoms:    
            atom1_acceptors = [acc for acc in self.hba_atoms if acc.idx == atom1.idx]    
            atom2_acceptors = [acc for acc in self.hba_atoms if acc.idx == atom2.idx]    
            atom1_donors = [don for don in self.hbd_atoms if don.d.idx == atom1.idx]    
            atom2_donors = [don for don in self.hbd_atoms if don.d.idx == atom2.idx]  

            # 检测氢键 - 使用精细的 hbonds 函数  
            if atom1_acceptors and atom2_donors:  
                hbond_results = hbonds(atom1_acceptors, atom2_donors)  
                if hbond_results:  
                    interaction_types.append('hydrogen_bond')  
            elif atom2_acceptors and atom1_donors:  
                hbond_results = hbonds(atom2_acceptors, atom1_donors)  
                if hbond_results:  
                    interaction_types.append('hydrogen_bond')  

        # 2. 检测疏水相互作用 - 使用精细的 hydrophobic_interactions 函数  
        if self.hydrophobic_atoms:  
            atom1_hydrophobic = [h for h in self.hydrophobic_atoms if h.idx == atom1.idx]  
            atom2_hydrophobic = [h for h in self.hydrophobic_atoms if h.idx == atom2.idx]  

            if atom1_hydrophobic and atom2_hydrophobic:  
                hydrophobic_results = hydrophobic_interactions(atom1_hydrophobic, atom2_hydrophobic)  
                if hydrophobic_results:  
                    interaction_types.append('hydrophobic')  

        # 3. 检测卤键 - 使用精细的 halogen_bonds 函数  
        if self.hal_acceptors and self.hal_donors:  
            atom1_hal_acc = [acc for acc in self.hal_acceptors if acc.o.idx == atom1.idx]  
            atom2_hal_acc = [acc for acc in self.hal_acceptors if acc.o.idx == atom2.idx]  
            atom1_hal_don = [don for don in self.hal_donors if don.x.idx == atom1.idx]  
            atom2_hal_don = [don for don in self.hal_donors if don.x.idx == atom2.idx]  

            if atom1_hal_acc and atom2_hal_don:  
                halogen_results = halogen_bonds(atom1_hal_acc, atom2_hal_don)  
                if halogen_results:  
                    interaction_types.append('halogen_bond')  
            elif atom2_hal_acc and atom1_hal_don:  
                halogen_results = halogen_bonds(atom2_hal_acc, atom1_hal_don)  
                if halogen_results:  
                    interaction_types.append('halogen_bond')  

        # 4. 检测盐桥 - 使用精细的 saltbridge 函数  
        if self.charged_groups:  
            pos_charged = [g for g in self.charged_groups if g.type == 'positive']  
            neg_charged = [g for g in self.charged_groups if g.type == 'negative']  

            # 检查原子是否属于带电基团  
            atom1_pos = [g for g in pos_charged if any(a.idx == atom1.idx for a in g.atoms)]  
            atom1_neg = [g for g in neg_charged if any(a.idx == atom1.idx for a in g.atoms)]  
            atom2_pos = [g for g in pos_charged if any(a.idx == atom2.idx for a in g.atoms)]  
            atom2_neg = [g for g in neg_charged if any(a.idx == atom2.idx for a in g.atoms)]  

            if (atom1_pos and atom2_neg) or (atom1_neg and atom2_pos):  
                # 使用相关的带电基团进行盐桥检测  
                relevant_pos = atom1_pos if atom1_pos else atom2_pos  
                relevant_neg = atom1_neg if atom1_neg else atom2_neg  
                saltbridge_results = saltbridge(relevant_pos, relevant_neg)  
                if saltbridge_results:  
                    interaction_types.append('salt_bridge')  

        return {  
            'interaction_types': interaction_types,  
            'distance': distance  
        } if interaction_types else None