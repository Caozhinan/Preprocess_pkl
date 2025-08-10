"""  
水桥和金属配位检测脚本  
基于PLIP (Protein-Ligand Interaction Profiler) 的核心算法  
集成到affincraft预处理流程中  
"""  
  
import itertools  
import math  
from collections import namedtuple, defaultdict  
from typing import List, Tuple, Dict, Any  
import numpy as np  
  
# 配置参数  
class WaterMetalConfig:  
    # 水桥检测参数  
    WATER_BRIDGE_MINDIST = 2.5  # 水分子氧原子与极性原子的最小距离  
    WATER_BRIDGE_MAXDIST = 4.1  # 水分子氧原子与极性原子的最大距离  
    WATER_BRIDGE_OMEGA_MIN = 71  # 受体-水氧-供体氢之间的最小角度  
    WATER_BRIDGE_OMEGA_MAX = 140  # 受体-水氧-供体氢之间的最大角度  
    WATER_BRIDGE_THETA_MIN = 100  # 水氧-供体氢-供体原子之间的最小角度  
      
    # 金属配位参数  
    METAL_DIST_MAX = 3.0  # 金属离子与配位原子的最大距离  
      
    # 通用参数  
    MIN_DIST = 0.5  # 所有距离阈值的最小距离  
  
config = WaterMetalConfig()  
  
# 几何计算函数  
def euclidean3d(coord1: Tuple[float, float, float], coord2: Tuple[float, float, float]) -> float:  
    """计算两点间的欧几里得距离"""  
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))  
  
def vector(coord1: Tuple[float, float, float], coord2: Tuple[float, float, float]) -> Tuple[float, float, float]:  
    """计算从coord1到coord2的向量"""  
    return tuple(b - a for a, b in zip(coord1, coord2))  
  
def vecangle(vec1: Tuple[float, float, float], vec2: Tuple[float, float, float]) -> float:  
    """计算两个向量之间的角度（度）"""  
    def dot_product(v1, v2):  
        return sum(a * b for a, b in zip(v1, v2))  
      
    def magnitude(v):  
        return math.sqrt(sum(x ** 2 for x in v))  
      
    mag1, mag2 = magnitude(vec1), magnitude(vec2)  
    if mag1 == 0 or mag2 == 0:  
        return 0.0  
      
    cos_angle = dot_product(vec1, vec2) / (mag1 * mag2)  
    cos_angle = max(-1.0, min(1.0, cos_angle))  # 防止数值误差  
    return math.degrees(math.acos(cos_angle))  
  
# 适配到现有系统的检测函数  
def detect_water_bridges_from_atoms(src_atom, tgt_atom, distance, molecule):  
    """  
    检测两个原子之间是否存在水桥相互作用  
    适配到classify_protein_spatial_edges的调用方式  
      
    Args:  
        src_atom: 源原子 (AtomInfo对象)  
        tgt_atom: 目标原子 (AtomInfo对象)  
        distance: 原子间距离  
        molecule: OpenBabel分子对象  
      
    Returns:  
        bool: 是否存在水桥相互作用  
    """  
    try:  
        # 检查距离是否在合理范围内  
        if distance > config.WATER_BRIDGE_MAXDIST * 2:  # 水桥可能涉及更长的距离  
            return False  
          
        # 简化的水桥检测逻辑  
        # 1. 检查是否有水分子参与  
        src_is_water = _is_water_atom(src_atom)  
        tgt_is_water = _is_water_atom(tgt_atom)  
          
        if src_is_water or tgt_is_water:  
            # 如果其中一个是水原子，检查另一个是否是极性原子  
            non_water_atom = tgt_atom if src_is_water else src_atom  
            if _is_polar_atom(non_water_atom):  
                return True  
          
        # 2. 检查是否可能通过水分子桥接  
        # 这需要更复杂的分子环境分析，这里提供简化版本  
        if _is_polar_atom(src_atom) and _is_polar_atom(tgt_atom):  
            # 如果两个都是极性原子且距离适中，可能存在水桥  
            if config.WATER_BRIDGE_MINDIST * 2 <= distance <= config.WATER_BRIDGE_MAXDIST * 2:  
                # 进一步检查分子中是否有水分子  
                if _has_nearby_water(src_atom, tgt_atom, molecule):  
                    return True  
          
        return False  
          
    except Exception as e:  
        # 如果检测失败，返回False  
        return False  
  
def detect_metal_complex_from_atoms(src_atom, tgt_atom, distance, molecule):  
    """  
    检测两个原子之间是否存在金属配位相互作用  
    适配到classify_protein_spatial_edges的调用方式  
      
    Args:  
        src_atom: 源原子 (AtomInfo对象)  
        tgt_atom: 目标原子 (AtomInfo对象)  
        distance: 原子间距离  
        molecule: OpenBabel分子对象  
      
    Returns:  
        bool: 是否存在金属配位相互作用  
    """  
    try:  
        # 检查距离是否在金属配位范围内  
        if distance > config.METAL_DIST_MAX:  
            return False  
          
        # 金属原子类型列表  
        metal_types = ['ZN', 'MG', 'CA', 'FE', 'MN', 'CU', 'NI', 'CO', 'CD', 'HG']  
          
        # 检查是否有金属原子参与  
        src_is_metal = _is_metal_atom(src_atom, metal_types)  
        tgt_is_metal = _is_metal_atom(tgt_atom, metal_types)  
          
        if src_is_metal or tgt_is_metal:  
            # 如果其中一个是金属原子，检查另一个是否是配位原子  
            non_metal_atom = tgt_atom if src_is_metal else src_atom  
            if _is_coordination_atom(non_metal_atom):  
                return True  
          
        return False  
          
    except Exception as e:  
        # 如果检测失败，返回False  
        return False  
  
# 辅助函数  
def _is_water_atom(atom):  
    """检查原子是否是水分子中的原子"""  
    try:  
        # 检查残基类型  
        if hasattr(atom, 'restype') and atom.restype in ['HOH', 'WAT', 'H2O']:  
            return True  
          
        # 检查原子类型和环境  
        if hasattr(atom, 'type') and atom.type == 'O':  
            # 进一步检查是否在水分子环境中  
            # 这里可以添加更复杂的逻辑  
            return True  
              
        return False  
    except:  
        return False  
  
def _is_polar_atom(atom):  
    """检查原子是否是极性原子（可参与氢键）"""  
    try:  
        if hasattr(atom, 'type'):  
            atom_type = atom.type.upper()  
            # 常见的极性原子类型  
            polar_types = ['O', 'N', 'S', 'F']  
            return atom_type in polar_types  
        return False  
    except:  
        return False  
  
def _is_metal_atom(atom, metal_types):  
    """检查原子是否是金属原子"""  
    try:  
        if hasattr(atom, 'type'):  
            atom_type = atom.type.upper()  
            return atom_type in metal_types  
        return False  
    except:  
        return False  
  
def _is_coordination_atom(atom):  
    """检查原子是否可以作为配位原子"""  
    try:  
        if hasattr(atom, 'type'):  
            atom_type = atom.type.upper()  
            # 常见的配位原子类型  
            coordination_types = ['O', 'N', 'S', 'CL', 'BR', 'I', 'F']  
            return atom_type in coordination_types  
        return False  
    except:  
        return False  
  
def _has_nearby_water(src_atom, tgt_atom, molecule):  
    """检查两个原子附近是否有水分子"""  
    try:  
        # 简化实现：检查分子中是否有水原子  
        for atom in molecule.atoms:  
            if _is_water_atom_simple(atom):  
                # 计算水原子到两个目标原子的距离  
                water_coord = np.array(atom.coords)  
                src_coord = np.array(src_atom.coords)  
                tgt_coord = np.array(tgt_atom.coords)  
                  
                dist_to_src = np.linalg.norm(water_coord - src_coord)  
                dist_to_tgt = np.linalg.norm(water_coord - tgt_coord)  
                  
                if (dist_to_src <= config.WATER_BRIDGE_MAXDIST and   
                    dist_to_tgt <= config.WATER_BRIDGE_MAXDIST):  
                    return True  
          
        return False  
    except:  
        return False  
  
def _is_water_atom_simple(atom):  
    """简化的水原子检测"""  
    try:  
        # 检查原子是否在水分子残基中  
        residue = atom.OBAtom.GetResidue()  
        if residue:  
            res_name = residue.GetName()  
            if res_name in ['HOH', 'WAT', 'H2O']:  
                return True  
        return False  
    except:  
        return False  