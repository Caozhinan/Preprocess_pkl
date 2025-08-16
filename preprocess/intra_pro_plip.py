#!/usr/bin/env python3  
"""  
基于PLIP和OpenBabel的非共价相互作用类型识别脚本  
用于分析蛋白质内部原子对的相互作用类型  
"""  
  
import itertools  
import numpy as np  
from collections import namedtuple  
from openbabel import pybel  
from openbabel.openbabel import OBAtomAtomIter  
  
# 配置参数 (基于PLIP的默认阈值)  
class Config:  
    MIN_DIST = 0.5  
    HYDROPH_DIST_MAX = 4.0  
    HBOND_DIST_MAX = 4.0  
    HBOND_DON_ANGLE_MIN = 100  
    PISTACK_DIST_MAX = 5.5  
    PISTACK_ANG_DEV = 30  
    PISTACK_OFFSET_MAX = 2.0  
    PICATION_DIST_MAX = 6.0  
    SALTBRIDGE_DIST_MAX = 5.5  
    HALOGEN_DIST_MAX = 4.0  
    HALOGEN_ACC_ANGLE = 120  
    HALOGEN_DON_ANGLE = 165  
    HALOGEN_ANGLE_DEV = 30  
    AROMATIC_PLANARITY = 5.0  
  
config = Config()  
  
def euclidean3d(coord1, coord2):  
    """计算两点间的欧几里得距离"""  
    return np.sqrt(sum([(coord1[i] - coord2[i])**2 for i in range(3)]))  
  
def vector(coord1, coord2):  
    """计算从coord1到coord2的向量"""  
    return np.array([coord2[i] - coord1[i] for i in range(3)])  
  
def vecangle(vec1, vec2):  
    """计算两个向量间的夹角(度)"""  
    cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))  
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  
    return np.degrees(np.arccos(cos_angle))  
  
def normalize_vector(v):  
    """向量归一化"""  
    norm = np.linalg.norm(v)  
    if norm == 0:  
        return v  
    return v / norm  
  
def projection(normal, center, point):  
    """将点投影到由法向量和中心定义的平面上"""  
    v = np.array(point) - np.array(center)  
    proj_length = np.dot(v, normal) / np.dot(normal, normal)  
    return np.array(point) - proj_length * np.array(normal)  
  
def centroid(coords):  
    """计算坐标的质心"""  
    return np.mean(coords, axis=0)  
  
def ring_is_planar(ring, atoms):  
    """检查环是否平面，基于PLIP的实现""" [1] 
    if len(atoms) < 4:  
        return True  
      
    # 计算每个原子到其邻居的法向量  
    normals = []  
    for i in range(len(atoms)):  
        prev_atom = atoms[i-1]  
        curr_atom = atoms[i]  
        next_atom = atoms[(i+1) % len(atoms)]  
          
        v1 = vector(curr_atom.coords, prev_atom.coords)  
        v2 = vector(curr_atom.coords, next_atom.coords)  
        normal = np.cross(v1, v2)  
        if np.linalg.norm(normal) > 0:  
            normals.append(normalize_vector(normal))  
      
    if len(normals) < 2:  
        return True  
      
    # 检查所有法向量之间的角度  
    for i in range(len(normals)):  
        for j in range(i+1, len(normals)):  
            angle = vecangle(normals[i], normals[j])  
            if angle > config.AROMATIC_PLANARITY and (180 - angle) > config.AROMATIC_PLANARITY:  
                return False  
    return True  
  
class AtomInfo:  
    """原子信息类，包装OpenBabel原子"""  
    def __init__(self, pybel_atom, idx=None):  
        self.atom = pybel_atom  
        self.idx = idx if idx is not None else pybel_atom.idx  
        self.coords = pybel_atom.coords  
        self.atomicnum = pybel_atom.atomicnum  
        self.type = pybel_atom.type  
        self.OBAtom = pybel_atom.OBAtom  
  
class RingInfo:  
    """芳香环信息类"""  
    def __init__(self, atoms, center, normal, ring_type='aromatic'):  
        self.atoms = atoms  
        self.center = center  
        self.normal = normal  
        self.type = ring_type  
  
def find_hydrophobic_atoms(molecule):  
    """识别疏水原子 - 基于PLIP的hydrophobic_atoms方法""" [2]   
    hydrophobic_atoms = []  
      
    for atom in molecule.atoms:  
        # 只考虑碳原子  
        if atom.atomicnum == 6:  
            # 检查邻居原子是否只有碳和氢  
            neighbor_nums = set()  
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):  
                neighbor_nums.add(neighbor.GetAtomicNum())  
              
            # 如果邻居只有碳(6)和氢(1)，则为疏水原子  
            if neighbor_nums.issubset({1, 6}):  
                hydrophobic_atoms.append(AtomInfo(atom))  
      
    return hydrophobic_atoms  
  
def find_hba(molecule):  
    """识别氢键受体 - 基于PLIP的find_hba方法""" [3] 
    acceptors = []  
      
    for atom in molecule.atoms:  
        # 排除卤素原子  
        if atom.atomicnum not in [9, 17, 35, 53]:  
            if atom.OBAtom.IsHbondAcceptor():  
                acceptors.append(AtomInfo(atom))  
      
    return acceptors  
  
def find_hbd(molecule, hydrophobic_atoms=None):  
    """识别氢键供体 - 基于PLIP的find_hbd方法""" [4]   
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
  
def find_rings(molecule):  
    """识别芳香环 - 基于PLIP的find_rings方法""" [5]   
    rings = []  
    aromatic_amino = ['TYR', 'TRP', 'HIS', 'PHE']  
      
    ring_candidates = molecule.OBMol.GetSSSR()  
      
    for ring in ring_candidates:  
        ring_atoms = []  
        for atom in molecule.atoms:  
            if ring.IsMember(atom.OBAtom):  
                ring_atoms.append(AtomInfo(atom))  
          
        ring_atoms = sorted(ring_atoms, key=lambda x: x.idx)  
          
        if 4 < len(ring_atoms) <= 6:  
            # 检查是否为芳香环  
            if ring.IsAromatic() or ring_is_planar(ring, ring_atoms):  
                # 计算环中心和法向量  
                ring_coords = [atom.coords for atom in ring_atoms]  
                center = centroid(ring_coords)  
                  
                # 使用三个原子计算法向量  
                if len(ring_atoms) >= 3:  
                    v1 = vector(ring_coords[0], ring_coords[2])  
                    v2 = vector(ring_coords[4 if len(ring_atoms) > 4 else 1], ring_coords[0])  
                    normal = normalize_vector(np.cross(v1, v2))  
                      
                    ring_info = RingInfo(ring_atoms, center, normal, f'{len(ring_atoms)}-membered')  
                    rings.append(ring_info)  
      
    return rings  
  
def find_charged_groups(molecule):  
    """识别带电基团 - 基于PLIP的find_charged方法""" [6]   
    data = namedtuple('charge', 'atoms type center')  
    charged_groups = []  
      
    # 遍历所有残基  
    for res in pybel.ob.OBResidueIter(molecule.OBMol):  
        contributing_atoms = []  
          
        # 正电荷：精氨酸、组氨酸、赖氨酸  
        if res.GetName() in ('ARG', 'HIS', 'LYS'):  
            for atom in pybel.ob.OBResidueAtomIter(res):  
                if atom.GetType().startswith('N') and res.GetAtomProperty(atom, 8):  # 侧链原子  
                    contributing_atoms.append(pybel.Atom(atom))  
              
            if contributing_atoms:  
                center = centroid([atom.coords for atom in contributing_atoms])  
                charged_groups.append(data(atoms=contributing_atoms, type='positive', center=center))  
          
        # 负电荷：天冬氨酸、谷氨酸  
        elif res.GetName() in ('GLU', 'ASP'):  
            for atom in pybel.ob.OBResidueAtomIter(res):  
                if atom.GetType().startswith('O') and res.GetAtomProperty(atom, 8):  # 侧链原子  
                    contributing_atoms.append(pybel.Atom(atom))  
              
            if contributing_atoms:  
                center = centroid([atom.coords for atom in contributing_atoms])  
                charged_groups.append(data(atoms=contributing_atoms, type='negative', center=center))  
      
    return charged_groups  
  
def find_halogen_acceptors(molecule):  
    """识别卤键受体 - 基于PLIP的find_hal方法""" [7] 
    data = namedtuple('hal_acceptor', 'o y')  
    acceptors = []  
      
    for atom in molecule.atoms:  
        # 氧、氮、硫原子  
        if atom.atomicnum in [8, 7, 16]:  
            neighbors = []  
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):  
                if neighbor.GetAtomicNum() in [6, 7, 15, 16]:  # C, N, P, S  
                    neighbors.append(pybel.Atom(neighbor))  
              
            if len(neighbors) == 1:  # 只有一个邻接原子  
                acceptors.append(data(o=AtomInfo(atom), y=neighbors[0]))  
      
    return acceptors  
  
def find_halogen_donors(molecule):  
    """识别卤键供体"""  
    data = namedtuple('hal_donor', 'x c')  
    donors = []  
      
    for atom in molecule.atoms:  
        # 卤素原子  
        if atom.atomicnum in [9, 17, 35, 53]:  # F, Cl, Br, I  
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):  
                if neighbor.GetAtomicNum() == 6:  # 连接到碳原子  
                    donors.append(data(x=AtomInfo(atom), c=pybel.Atom(neighbor)))  
                    break  
      
    return donors  
  
# 主要检测函数的完整实现  
  
def hydrophobic_interactions(atom_set_a, atom_set_b):  
    """检测疏水相互作用 - 基于PLIP的hydrophobic_interactions函数""" [8]   
    data = namedtuple('hydroph_interaction', 'atom1 atom2 distance type')  
    pairings = []  
      
    for a, b in itertools.product(atom_set_a, atom_set_b):  
        if a.idx == b.idx:  
            continue  
        distance = euclidean3d(a.coords, b.coords)  
        if config.MIN_DIST < distance < config.HYDROPH_DIST_MAX:  
            contact = data(atom1=a, atom2=b, distance=distance, type='hydrophobic')  
            pairings.append(contact)  
      
    return pairings  


def hbonds(acceptors, donors):  
    """检测氢键 - 基于PLIP的hbonds函数"""  
    data = namedtuple('hbond', 'acceptor donor hydrogen distance_ah distance_ad angle type')  
    pairings = []  
      
    for acc, don in itertools.product(acceptors, donors):  
        # 计算供体-受体距离  
        dist_ad = euclidean3d(acc.coords, don.d.coords)  
        if not config.MIN_DIST < dist_ad < config.HBOND_DIST_MAX:  
            continue  
          
        # 计算氢-受体距离  
        dist_ah = euclidean3d(acc.coords, don.h.coords)  
          
        # 计算供体角度 D-H...A  
        vec1 = vector(don.h.coords, don.d.coords)  
        vec2 = vector(don.h.coords, acc.coords)  
        angle = vecangle(vec1, vec2)  
          
        if angle > config.HBOND_DON_ANGLE_MIN:  
            contact = data(acceptor=acc, donor=don.d, hydrogen=don.h,  
                         distance_ah=dist_ah, distance_ad=dist_ad,   
                         angle=angle, type='hydrogen_bond')  
            pairings.append(contact)  
      
    return pairings

def pistacking(rings_a, rings_b):  
    """检测π-π堆积"""  
    data = namedtuple('pistacking', 'ring1 ring2 distance angle offset type')  
    pairings = []  
      
    for ring1, ring2 in itertools.product(rings_a, rings_b):  
        distance = euclidean3d(ring1.center, ring2.center)  
        if config.MIN_DIST < distance < config.PISTACK_DIST_MAX:  
            # 计算环平面间的角度  
            angle = vecangle(ring1.normal, ring2.normal)  
            angle = min(angle, 180 - angle)  # 取较小角度  
              
            # 计算偏移  
            proj1 = projection(ring2.normal, ring2.center, ring1.center)  
            proj2 = projection(ring1.normal, ring1.center, ring2.center)  
            offset = min(euclidean3d(proj1, ring2.center), euclidean3d(proj2, ring1.center))  
              
            # 判断堆积类型  
            if angle < config.PISTACK_ANG_DEV and offset < config.PISTACK_OFFSET_MAX:  
                stack_type = 'parallel'  
                contact = data(ring1=ring1, ring2=ring2, distance=distance,  
                             angle=angle, offset=offset, type=f'pi_stacking_{stack_type}')  
                pairings.append(contact)  
            elif abs(90 - angle) < config.PISTACK_ANG_DEV and offset < config.PISTACK_OFFSET_MAX:  
                stack_type = 'perpendicular'  
                contact = data(ring1=ring1, ring2=ring2, distance=distance,  
                             angle=angle, offset=offset, type=f'pi_stacking_{stack_type}')  
                pairings.append(contact)  
      
    return pairings  
  
def saltbridge(pos_charged, neg_charged):  
    """检测盐桥"""  
    data = namedtuple('saltbridge', 'positive negative distance type')  
    pairings = []  
      
    for pos, neg in itertools.product(pos_charged, neg_charged):  
        distance = euclidean3d(pos.center, neg.center)  
        if config.MIN_DIST < distance < config.SALTBRIDGE_DIST_MAX:  
            contact = data(positive=pos, negative=neg, distance=distance, type='salt_bridge')  
            pairings.append(contact)  
      
    return pairings  
  
def pication(rings, charged_groups):  
    """检测π-阳离子相互作用"""  
    data = namedtuple('pication', 'ring charge distance offset type')  
    pairings = []  
      
    for ring, charge in itertools.product(rings, charged_groups):  
        distance = euclidean3d(ring.center, charge.center)  
        if config.MIN_DIST < distance < config.PICATION_DIST_MAX:  
            # 计算偏移  
            proj = projection(ring.normal, ring.center, charge.center)  
            offset = euclidean3d(proj, ring.center)  
              
            if offset < config.PISTACK_OFFSET_MAX:  
                contact = data(ring=ring, charge=charge, distance=distance,  
                             offset=offset, type='pi_cation')  
                pairings.append(contact)  
      
    return pairings  
  
def halogen_bonds(acceptors, donors):  
    """检测卤键"""  
    data = namedtuple('halogen_bond', 'acceptor donor distance don_angle acc_angle type')  
    pairings = []  
      
    for acc, don in itertools.product(acceptors, donors):  
        distance = euclidean3d(acc.o.coords, don.x.coords)  
        if config.MIN_DIST < distance < config.HALOGEN_DIST_MAX:  
            # 计算角度  
            vec1 = vector(acc.o.coords, acc.y.coords)  
            vec2 = vector(acc.o.coords, don.x.coords)  
            acc_angle = vecangle(vec1, vec2)  
              
            vec3 = vector(don.x.coords, acc.o.coords)  
            vec4 = vector(don.x.coords, don.c.coords)  
            don_angle = vecangle(vec3, vec4)  
              
            # 检查角度条件  
            if (config.HALOGEN_ACC_ANGLE - config.HALOGEN_ANGLE_DEV < acc_angle <   
                config.HALOGEN_ACC_ANGLE + config.HALOGEN_ANGLE_DEV and  
                config.HALOGEN_DON_ANGLE - config.HALOGEN_ANGLE_DEV < don_angle <   
                config.HALOGEN_DON_ANGLE + config.HALOGEN_ANGLE_DEV):  
                  
                contact = data(acceptor=acc, donor=don, distance=distance,  
                             don_angle=don_angle, acc_angle=acc_angle, type='halogen_bond')  
                pairings.append(contact)  
      
    return pairings

class InteractionAnalyzer:  
    """相互作用分析器主类 - 使用精细的PLIP检测函数"""  
      
    def __init__(self):  
        self.interactions = []  
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
      
    def analyze_atom_pair(self, atom1, atom2, distance=None):  
        """分析单个原子对的相互作用类型 - 使用精细检测函数"""  
        if distance is None:  
            distance = euclidean3d(atom1.coords, atom2.coords)  
          
        if distance > 5.0:  # 超过5Å的不分析  
            return None  
          
        interaction_types = []  
          
        # 1. 检测氢键 - 使用精细的hbonds函数  
        if self.hba_atoms and self.hbd_atoms:  
            atom1_acceptors = [acc for acc in self.hba_atoms if acc.idx == atom1.idx]  
            atom2_acceptors = [acc for acc in self.hba_atoms if acc.idx == atom2.idx]  
            atom1_donors = [don for don in self.hbd_atoms if don.d.idx == atom1.idx]  
            atom2_donors = [don for don in self.hbd_atoms if don.d.idx == atom2.idx]  
              
            # 检测氢键  
            if atom1_acceptors and atom2_donors:  
                hbond_results = hbonds(atom1_acceptors, atom2_donors)  
                if hbond_results:  
                    interaction_types.append('hydrogen_bond')  
            elif atom2_acceptors and atom1_donors:  
                hbond_results = hbonds(atom2_acceptors, atom1_donors)  
                if hbond_results:  
                    interaction_types.append('hydrogen_bond')  
          
        # 2. 检测疏水相互作用 - 使用精细的hydrophobic_interactions函数  
        if self.hydrophobic_atoms:  
            atom1_hydrophobic = [h for h in self.hydrophobic_atoms if h.idx == atom1.idx]  
            atom2_hydrophobic = [h for h in self.hydrophobic_atoms if h.idx == atom2.idx]  
              
            if atom1_hydrophobic and atom2_hydrophobic:  
                hydrophobic_results = hydrophobic_interactions(atom1_hydrophobic, atom2_hydrophobic)  
                if hydrophobic_results:  
                    interaction_types.append('hydrophobic')  
          
        # 3. 检测卤键 - 使用精细的halogen_bonds函数  
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
          
        # 4. 检测盐桥 - 使用精细的saltbridge函数  
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
            'atom1': atom1,  
            'atom2': atom2,  
            'distance': distance,  
            'interaction_types': interaction_types  
        }  
      
    def analyze_molecule_interactions(self, molecule):  
        """分析整个分子的相互作用（包括π-π堆积和π-阳离子）- 使用精细检测函数"""  
        if not self.molecule:  
            self.set_molecule(molecule)  
          
        results = []  
          
        # 检测π-π堆积 - 使用精细的pistacking函数  
        if self.rings:  
            pi_stacks = pistacking(self.rings, self.rings)  
            for stack in pi_stacks:  
                results.append({  
                    'interaction_type': 'pi_stacking',  
                    'ring1': stack.ring1,  
                    'ring2': stack.ring2,  
                    'distance': stack.distance,  
                    'type': stack.type  
                })  
          
        # 检测π-阳离子相互作用 - 使用精细的pication函数  
        if self.rings and self.charged_groups:  
            pos_charged = [g for g in self.charged_groups if g.type == 'positive']  
            pi_cations = pication(self.rings, pos_charged)  
            for pication_int in pi_cations:  
                results.append({  
                    'interaction_type': 'pi_cation',  
                    'ring': pication_int.ring,  
                    'charge': pication_int.charge,  
                    'distance': pication_int.distance  
                })  
          
        return results
  
# 使用示例  
if __name__ == "__main__":  
    analyzer = InteractionAnalyzer()  
    print("完整的非共价相互作用分析脚本已准备就绪")  
    print("支持检测：疏水相互作用、氢键、盐桥、卤键、π-π堆积、π-阳离子相互作用")