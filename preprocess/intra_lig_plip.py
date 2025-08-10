#!/usr/bin/env python3  
"""  
配体内部相互作用分析模块  
直接复用 intra_pro_plip.py 的所有类和函数定义，适配配体分子  
"""  
  
import numpy as np  
from openbabel import pybel  
from typing import List, Dict, Tuple, Optional, Set  
import math  
  
# 直接复用 intra_pro_plip.py 中的所有类和函数定义  
class AtomInfo:  
    """原子信息类，存储原子的基本属性"""  
    def __init__(self, atom):  
        self.atom = atom  
        self.idx = atom.idx  
        self.atomicnum = atom.atomicnum  
        self.coords = np.array([atom.coords[0], atom.coords[1], atom.coords[2]])  
        self.element = atom.OBAtom.GetAtomicNum()  
        self.formal_charge = atom.formalcharge  
        self.partial_charge = getattr(atom.OBAtom, 'GetPartialCharge', lambda: 0.0)()  
  
class Config:  
    """配置参数类 - 与 intra_pro_plip.py 保持一致"""  
    HYDROPHOBIC_DIST_MAX = 4.0  
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
    WATER_BRIDGE_MINDIST = 2.5  
    WATER_BRIDGE_MAXDIST = 4.0  
    METAL_DIST_MAX = 3.0  
  
# 直接复用 intra_pro_plip.py 中的所有功能函数  
def find_hydrophobic_atoms(molecule) -> Set[int]:  
    """识别疏水原子 - 复用蛋白质分析逻辑"""  
    hydrophobic_atoms = set()  
      
    for atom in molecule.atoms:  
        if atom.atomicnum == 6:  # Carbon  
            polar_neighbors = 0  
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):  
                if neighbor.GetAtomicNum() in [7, 8, 9, 15, 16, 17, 35, 53]:  
                    polar_neighbors += 1  
              
            if polar_neighbors <= 1:  
                hydrophobic_atoms.add(atom.idx)  
        elif atom.atomicnum == 16:  # Sulfur  
            hydrophobic_atoms.add(atom.idx)  
      
    return hydrophobic_atoms  
  
def find_hba(molecule) -> Set[int]:  
    """识别氢键受体原子 - 复用蛋白质分析逻辑"""  
    hba_atoms = set()  
      
    for atom in molecule.atoms:  
        if atom.atomicnum in [7, 8, 9]:  # N, O, F  
            if atom.atomicnum == 7:  # Nitrogen  
                valence = len(list(pybel.ob.OBAtomAtomIter(atom.OBAtom)))  
                if valence < 4:  
                    hba_atoms.add(atom.idx)  
            elif atom.atomicnum == 8:  # Oxygen  
                hba_atoms.add(atom.idx)  
            elif atom.atomicnum == 9:  # Fluorine  
                hba_atoms.add(atom.idx)  
      
    return hba_atoms  
  
def find_hbd(molecule, hydrophobic_atoms: Set[int]) -> Set[int]:  
    """识别氢键供体原子 - 复用蛋白质分析逻辑"""  
    hbd_atoms = set()  
      
    for atom in molecule.atoms:  
        if atom.atomicnum == 1:  # Hydrogen  
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):  
                if neighbor.GetAtomicNum() in [7, 8]:  # N, O  
                    hbd_atoms.add(atom.idx)  
                    break  
      
    return hbd_atoms  
  
class RingInfo:  
    """芳香环信息类 - 与 intra_pro_plip.py 保持一致"""  
    def __init__(self, atoms, center, normal, ring_type='aromatic'):  
        self.atoms = atoms  # 这里应该是 AtomInfo 对象列表  
        self.center = center  
        self.normal = normal  
        self.type = ring_type  
  
def find_rings(molecule):  
    """识别芳香环 - 返回 RingInfo 对象列表"""  
    rings = []  
      
    ring_candidates = molecule.OBMol.GetSSSR()  
      
    for ring in ring_candidates:  
        ring_atoms = []  
        for atom in molecule.atoms:  
            if ring.IsMember(atom.OBAtom):  
                ring_atoms.append(AtomInfo(atom))  # 创建 AtomInfo 对象  
          
        ring_atoms = sorted(ring_atoms, key=lambda x: x.idx)  
          
        if 4 < len(ring_atoms) <= 6:  
            if ring.IsAromatic():  
                ring_coords = [atom.coords for atom in ring_atoms]  
                center = np.mean(ring_coords, axis=0)  
                  
                if len(ring_atoms) >= 3:  
                    v1 = np.array(ring_coords[0]) - np.array(ring_coords[2])  
                    v2 = np.array(ring_coords[4 if len(ring_atoms) > 4 else 1]) - np.array(ring_coords[0])  
                    normal = np.cross(v1, v2)  
                    if np.linalg.norm(normal) > 0:  
                        normal = normal / np.linalg.norm(normal)  
                    else:  
                        normal = np.array([0, 0, 1])  
                  
                # 创建 RingInfo 对象  
                ring_obj = RingInfo(ring_atoms, center, normal)  
                rings.append(ring_obj)  
      
    return rings  
  
def find_charged_groups(molecule) -> Dict[str, List[int]]:  
    """识别带电基团 - 复用蛋白质分析逻辑"""  
    charged_groups = {'positive': [], 'negative': []}  
      
    for atom in molecule.atoms:  
        if atom.formalcharge > 0:  
            charged_groups['positive'].append(atom.idx)  
        elif atom.formalcharge < 0:  
            charged_groups['negative'].append(atom.idx)  
          
        if atom.atomicnum == 7:  # Nitrogen  
            neighbors = list(pybel.ob.OBAtomAtomIter(atom.OBAtom))  
            if len(neighbors) == 4:  
                charged_groups['positive'].append(atom.idx)  
          
        elif atom.atomicnum == 8:  # Oxygen  
            neighbors = list(pybel.ob.OBAtomAtomIter(atom.OBAtom))  
            for neighbor in neighbors:  
                if neighbor.GetAtomicNum() == 6:  
                    carbon_neighbors = list(pybel.ob.OBAtomAtomIter(neighbor))  
                    oxygen_count = sum(1 for n in carbon_neighbors if n.GetAtomicNum() == 8)  
                    if oxygen_count >= 2:  
                        charged_groups['negative'].append(atom.idx)  
                        break  
      
    return charged_groups  
  
def find_halogen_acceptors(molecule) -> Set[int]:  
    """识别卤键受体原子 - 复用蛋白质分析逻辑"""  
    halogen_acceptors = set()  
      
    for atom in molecule.atoms:  
        if atom.atomicnum in [7, 8, 16]:  # N, O, S  
            halogen_acceptors.add(atom.idx)  
      
    return halogen_acceptors  
  
def find_halogen_donors(molecule) -> Set[int]:  
    """识别卤键供体原子 - 复用蛋白质分析逻辑"""  
    halogen_donors = set()  
      
    for atom in molecule.atoms:  
        if atom.atomicnum in [9, 17, 35, 53]:  # F, Cl, Br, I  
            halogen_donors.add(atom.idx)  
      
    return halogen_donors  
  
def find_closest_atom_by_coord(coord: np.ndarray, molecule) -> Optional:  
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
    """相互作用分析器 - 复用蛋白质分析逻辑"""  
      
    def __init__(self):  
        self.config = Config()  
      
    def analyze_atom_pair(self, atom1, atom2, distance: float) -> Optional[Dict]:  
        """分析原子对之间的相互作用"""  
        interaction_types = []  
          
        if self._is_hydrogen_bond(atom1, atom2, distance):  
            interaction_types.append('hydrogen_bond')  
          
        if self._is_hydrophobic_interaction(atom1, atom2, distance):  
            interaction_types.append('hydrophobic')  
          
        if self._is_salt_bridge(atom1, atom2, distance):  
            interaction_types.append('salt_bridge')  
          
        if self._is_halogen_bond(atom1, atom2, distance):  
            interaction_types.append('halogen_bond')  
          
        if interaction_types:  
            return {  
                'interaction_types': interaction_types,  
                'distance': distance  
            }  
          
        return None  
      
    def _is_hydrogen_bond(self, atom1, atom2, distance: float) -> bool:  
        """检测氢键"""  
        if distance > self.config.HBOND_DIST_MAX:  
            return False  
          
        h_atom = None  
        heavy_atom = None  
          
        if atom1.atomicnum == 1:  
            h_atom = atom1  
            heavy_atom = atom2  
        elif atom2.atomicnum == 1:  
            h_atom = atom2  
            heavy_atom = atom1  
        else:  
            return False  
          
        if heavy_atom.atomicnum not in [7, 8, 9]:  
            return False  
          
        for neighbor in pybel.ob.OBAtomAtomIter(h_atom.OBAtom):  
            if neighbor.GetAtomicNum() in [7, 8]:  
                return True  
          
        return False  
      
    def _is_hydrophobic_interaction(self, atom1, atom2, distance: float) -> bool:  
        """检测疏水相互作用"""  
        if distance > self.config.HYDROPHOBIC_DIST_MAX:  
            return False  
          
        hydrophobic_elements = [6, 16]  # C, S  
          
        return (atom1.atomicnum in hydrophobic_elements and   
                atom2.atomicnum in hydrophobic_elements)  
      
    def _is_salt_bridge(self, atom1, atom2, distance: float) -> bool:  
        """检测盐桥"""  
        if distance > self.config.SALTBRIDGE_DIST_MAX:  
            return False  
          
        charge1 = atom1.formalcharge  
        charge2 = atom2.formalcharge  
          
        return (charge1 * charge2 < 0)  
      
    def _is_halogen_bond(self, atom1, atom2, distance: float) -> bool:  
        """检测卤键"""  
        if distance > self.config.HALOGEN_DIST_MAX:  
            return False  
          
        halogen_elements = [9, 17, 35, 53]  # F, Cl, Br, I  
        acceptor_elements = [7, 8, 16]  # N, O, S  
          
        return ((atom1.atomicnum in halogen_elements and atom2.atomicnum in acceptor_elements) or  
                (atom2.atomicnum in halogen_elements and atom1.atomicnum in acceptor_elements))  
  
def check_ring_interactions(atom1, atom2, distance: float, rings: List,   
                          charged_groups: Dict[str, List[int]]) -> str:  
    """检查环相关的相互作用 - 复用蛋白质分析逻辑"""  
      
    ring1 = None  
    ring2 = None  
      
    for ring in rings:  
        if atom1.idx in ring.atoms:  
            ring1 = ring  
        if atom2.idx in ring.atoms:  
            ring2 = ring  
      
    # π-π堆积检测  
    if ring1 and ring2 and ring1 != ring2:  
        if distance <= Config.PISTACK_DIST_MAX:  
            angle_between_normals = np.arccos(np.clip(np.dot(ring1.normal, ring2.normal), -1.0, 1.0))  
            angle_deg = np.degrees(angle_between_normals)  
              
            if angle_deg <= Config.PISTACK_ANG_DEV or abs(angle_deg - 90) <= Config.PISTACK_ANG_DEV:  
                center_distance = np.linalg.norm(ring1.center - ring2.center)  
                if center_distance <= Config.PISTACK_OFFSET_MAX:  
                    return 'pi_stacking'  
      
    # π-阳离子相互作用检测  
    if distance <= Config.PICATION_DIST_MAX:  
        atom1_in_ring = any(atom1.idx in ring.atoms for ring in rings)  
        atom2_in_ring = any(atom2.idx in ring.atoms for ring in rings)  
          
        atom1_charged = (atom1.idx in charged_groups['positive'] or   
                        atom1.idx in charged_groups['negative'])  
        atom2_charged = (atom2.idx in charged_groups['positive'] or   
                        atom2.idx in charged_groups['negative'])  
          
        if (atom1_in_ring and atom2_charged) or (atom2_in_ring and atom1_charged):  
            return 'pi_cation'  
      
    return 'others'

def detect_water_bridges_from_atoms(atom1, atom2, distance, molecule):  
    """  
    检测配体内部的水桥相互作用  
    直接调用现有的 water_metal_detection 模块  
    """  
    # 直接导入并调用现有的水桥检测函数  
    from water_metal_detection import detect_water_bridges_from_atoms as detect_water_original  
    return detect_water_original(atom1, atom2, distance, molecule)  
  
def detect_metal_complex_from_atoms(atom1, atom2, distance, molecule):  
    """  
    检测配体内部的金属配位相互作用  
    直接调用现有的 water_metal_detection 模块  
    """  
    # 直接导入并调用现有的金属配位检测函数  
    from water_metal_detection import detect_metal_complex_from_atoms as detect_metal_original  
    return detect_metal_original(atom1, atom2, distance, molecule)  
  
def classify_ligand_interactions(lig_sei, lig_sea, lig_coord, ligand_file):  
    """  
    为配体内部空间边分配相互作用类型  
      
    Args:  
        lig_sei: 配体空间边索引  
        lig_sea: 配体空间边属性  
        lig_coord: 配体原子坐标  
        ligand_file: 配体SDF文件路径  
      
    Returns:  
        tuple: (classified_edge_index, classified_edge_attr)  
    """  
    if len(lig_sei) == 0:  
        return lig_sei, lig_sea  
      
    try:  
        # 使用OpenBabel加载配体结构  
        molecule = pybel.readfile("sdf", ligand_file).__next__()  
          
        # 创建相互作用分析器  
        analyzer = InteractionAnalyzer()  
          
        # 识别各种功能基团  
        hydrophobic_atoms = find_hydrophobic_atoms(molecule)  
        hba_atoms = find_hba(molecule)  
        hbd_atoms = find_hbd(molecule, hydrophobic_atoms)  
        rings = find_rings(molecule)  
        charged_groups = find_charged_groups(molecule)  
        hal_acceptors = find_halogen_acceptors(molecule)  
        hal_donors = find_halogen_donors(molecule)  
          
        # 为每条空间边分配类型  
        classified_edge_attr = []  
          
        for i, (src_idx, tgt_idx) in enumerate(lig_sei.T):  
            src_coord = lig_coord[src_idx]  
            tgt_coord = lig_coord[tgt_idx]  
            distance = np.linalg.norm(src_coord - tgt_coord)  
              
            # 找到最接近的原子  
            src_atom = find_closest_atom_by_coord(src_coord, molecule)  
            tgt_atom = find_closest_atom_by_coord(tgt_coord, molecule)  
              
            if src_atom and tgt_atom:  
                # 使用InteractionAnalyzer分析原子对  
                interaction_result = analyzer.analyze_atom_pair(src_atom, tgt_atom, distance)  
                  
                if interaction_result and interaction_result['interaction_types']:  
                    # 取第一个检测到的相互作用类型  
                    interaction_type = interaction_result['interaction_types'][0]  
                      
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
                    # 检查π-π堆积和π-阳离子相互作用  
                    edge_type_name = check_ring_interactions(  
                        src_atom, tgt_atom, distance, rings, charged_groups  
                    )  
                      
                    # 检查水桥相互作用  
                    if edge_type_name == 'others' and detect_water_bridges_from_atoms(src_atom, tgt_atom, distance, molecule):  
                        edge_type_name = 'water_bridges'  
                      
                    # 检查金属配位相互作用  
                    if edge_type_name == 'others' and detect_metal_complex_from_atoms(src_atom, tgt_atom, distance, molecule):  
                        edge_type_name = 'metal_complexes'  
            else:  
                edge_type_name = 'others'  
              
            # 获取边类型编码  
            from preprocess import (HYDROGEN_BOND_EDGE, HYDROPHOBIC_EDGE, PI_STACKING_EDGE,   
                                  PI_CATION_EDGE, SALT_BRIDGE_EDGE, WATER_BRIDGE_EDGE,   
                                  HALOGEN_BOND_EDGE, METAL_COMPLEX_EDGE, OTHERS_EDGE)  
              
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
          
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)  
          
    except Exception as e:  
        print(f"配体相互作用分析失败: {e}")  
        # 如果分析失败，将所有边设置为OTHERS类型  
        from preprocess import OTHERS_EDGE  
        classified_edge_attr = []  
        for i, (src_idx, tgt_idx) in enumerate(lig_sei.T):  
            distance = np.linalg.norm(lig_coord[src_idx] - lig_coord[tgt_idx])  
            edge_feature = OTHERS_EDGE + [distance]  
            classified_edge_attr.append(edge_feature)  
        classified_edge_attr = np.array(classified_edge_attr, dtype=np.float32)  
      
    return lig_sei, classified_edge_attr