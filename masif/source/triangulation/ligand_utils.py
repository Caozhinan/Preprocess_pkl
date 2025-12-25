# 导入Python的shutil模块，用于文件操作（如复制文件）
import shutil
# 导入StringIO和BytesIO，用于字符串和字节流的内存操作
from io import StringIO, BytesIO
# 导入RDKit的化学分子处理库
from rdkit import Chem
from rdkit.Chem import AllChem , rdmolfiles
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
# 导入OpenBabel库用于分子格式转换
from openbabel import openbabel
# 导入ProDy用于PDB文件生物分子结构处理
import prody
# 设置ProDy的输出级别为最低（不显示日志信息）
import os  
from rdkit import Chem 
prody.confProDy(verbosity='none')


def neutralize_atoms(mol):
    # 定义SMARTS模式，匹配带正或负电荷的原子（排除某些情况）
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    # 在分子中查找所有符合SMARTS模式的原子
    at_matches = mol.GetSubstructMatches(pattern)
    # 提取所有匹配原子的索引
    at_matches_list = [y[0] for y in at_matches]
    # 如果有匹配，进行中和处理
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()    # 获取原子的形式电荷
            hcount = atom.GetTotalNumHs()   # 获取原子的氢原子数
            atom.SetFormalCharge(0)         # 电荷归零
            atom.SetNumExplicitHs(hcount - chg)  # 调整氢原子数
            atom.UpdatePropertyCache()      # 更新原子属性缓存
    # 返回一个新的分子对象（副本）
    mol_out = Chem.Mol(mol)
    return mol_out


def amide_to_single_bond(mol2_outfile):
    """
    修正mol2文件，将所有 am/ar 类型的键改为单键（1），并修正 .co2 原子类型
    """
    mol2_new = []
    bond_record = False
    atom_record = False
    # print("DEBUG: Using modified amide_to_single_bond function")

    try:
        with open(mol2_outfile, 'r') as f:
            for line in f:
                # 去除行尾换行
                line = line.rstrip('\n').rstrip('\r')

                # 区块切换
                if line.startswith('@<TRIPOS>'):
                    if line.startswith('@<TRIPOS>BOND'):
                        bond_record = True
                        atom_record = False
                    elif line.startswith('@<TRIPOS>ATOM'):
                        atom_record = True
                        bond_record = False
                    else:
                        bond_record = False
                        atom_record = False
                    mol2_new.append(line)
                    continue

                # 处理原子行
                if atom_record and line.strip() and not line.startswith('@<TRIPOS>'):
                    parts = line.split()
                    if len(parts) >= 6 and '.co2' in parts[5]:
                        parts[5] = parts[5].replace('.co2', '.2')
                        line = ' '.join(parts)
                    else:
                        line = ' '.join(parts)
                    mol2_new.append(line)
                    continue

                # 处理键行
                if bond_record and line.strip():
                    parts = line.split()
                    if len(parts) >= 4:
                        bond_id, origin_atom_id, target_atom_id, bond_type = parts[:4]
                        if bond_type.lower() in ('am', 'ar'):
                            bond_type = '1'
                        # 保留后面额外字段（如立体化学信息）
                        line = ' '.join([bond_id, origin_atom_id, target_atom_id, bond_type] + parts[4:])
                    else:
                        print(f"Warning: Skipping malformed bond line: '{line}'")
                    mol2_new.append(line)
                    continue

                # 其他行（注释、空行等）
                mol2_new.append(line)

        # 覆盖写回，注意用真换行
        with open(mol2_outfile, 'w') as f:
            f.write('\n'.join(mol2_new) + '\n')

    except Exception as e:
        print(f"Error processing MOL2 file {mol2_outfile}: {e}")
        raise

def _clean_mol2_format(mol2_file):
    """
    清理mol2格式：移除转义符、多余空白、非标准区块，仅保留正规Tripos区块
    """
    valid_blocks = {'@<TRIPOS>MOLECULE', '@<TRIPOS>ATOM', '@<TRIPOS>BOND', '@<TRIPOS>SUBSTRUCTURE'}
    output = []
    keep = False
    try:
        with open(mol2_file, 'r') as f:
            for line in f:
                line = line.replace('\\n', '').replace('\\t', '')
                line = line.rstrip('\n').rstrip('\r')
                if line.startswith('@<TRIPOS>'):
                    keep = line in valid_blocks
                    if keep:
                        output.append(line)
                    continue
                if keep:
                    if line.strip():
                        # 保证所有分隔符为单空格
                        split_line = line.split()
                        output.append(' '.join(split_line))
        # 写回
        with open(mol2_file, 'w') as f:
            f.write('\n'.join(output) + '\n')
    except Exception as e:
        print(f"Error cleaning MOL2 format: {e}")

def add_substructure_section(mol2_file):
    """
    补全/修正SUBSTRUCTURE区块（如果已经存在则不重复添加）
    """
    with open(mol2_file, 'r') as f:
        content = f.read()
    if '@<TRIPOS>SUBSTRUCTURE' not in content:
        with open(mol2_file, 'a') as f:
            f.write('@<TRIPOS>SUBSTRUCTURE\n')
            f.write('1 UNK 1 GROUP **** **** 0 ROOT\n')


def sdf_to_mol2(sdf_file, mol2_outfile):
    """
    Convert SDF file to MOL2 format using OpenBabel with format cleaning
    """
    try:
        from openbabel import openbabel

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "mol2")

        obmol = openbabel.OBMol()
        if not obConversion.ReadFile(obmol, sdf_file):
            print(f"Failed to read SDF file: {sdf_file}")
            return False

        if not obConversion.WriteFile(obmol, mol2_outfile):
            print(f"Failed to write MOL2 file: {mol2_outfile}")
            return False

        # 清理格式
        _clean_mol2_format(mol2_outfile)
        # 补全substructure
        add_substructure_section(mol2_outfile)
        # 修正键类型/原子类型
        amide_to_single_bond(mol2_outfile)
        # *** 保证原子名唯一 ***
        fix_atom_names(mol2_outfile)

        # print(f"[OK] SDF to MOL2 conversion and cleaning done: {mol2_outfile}")
        return True

    except Exception as e:
        print(f"Error during SDF to MOL2 conversion: {e}")
        return False

def fix_atom_names(mol2_path):
    """
    保证mol2文件ATOM区块每个原子名字唯一（如C1、C2、O1、H3...），避免pdb2pqr等工具报重复名。
    """
    new_lines = []
    atom_idx = {}
    in_atom = False
    with open(mol2_path, 'r') as f:
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom = True
                new_lines.append(line.rstrip())
                continue
            elif line.startswith('@<TRIPOS>'):
                in_atom = False
                new_lines.append(line.rstrip())
                continue
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    base = parts[1]
                    idx = atom_idx.get(base, 0) + 1
                    atom_idx[base] = idx
                    parts[1] = f"{base}{idx}"
                    line = ' '.join(parts)
            new_lines.append(line.rstrip())
    with open(mol2_path, 'w') as f:
        f.write('\n'.join(new_lines) + '\n')
    # print(f"[OK] Atom names fixed in {mol2_path}")
    
def extract_ligand(ligand_sdf_file="ligand.sdf"):    
    """    
    直接从SDF文件中提取配体分子    
    
    Args:    
        ligand_sdf_file: 配体SDF文件路径 (默认: "ligand.sdf")    
    
    Returns:    
        RDKit分子对象，如果失败则返回None    
    """    
    
    # 检查SDF文件是否存在    
    if not os.path.exists(ligand_sdf_file):    
        print(f"DEBUG: ERROR - SDF file does not exist: {ligand_sdf_file}")    
        return None    
    
    try:    
        # 修改：使用更宽松的参数读取分子  
        rdmol = Chem.MolFromMolFile(ligand_sdf_file, sanitize=False, removeHs=False)    
          
        if rdmol is not None:  
            # 手动进行基本的分子清理，但跳过价态检查  
            Chem.SanitizeMol(rdmol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)  
              
            # 为原子添加PDB残基信息  
            element_counts = {}    
            for i, atom in enumerate(rdmol.GetAtoms()):    
                info = Chem.AtomPDBResidueInfo()    
                element = atom.GetSymbol()    
                if element not in element_counts:    
                    element_counts[element] = 0    
                element_counts[element] += 1    
                if element_counts[element] == 1:    
                    atom_name = element    
                else:    
                    atom_name = f"{element}{element_counts[element]}"    
                info.SetName(atom_name)    
                info.SetResidueName("UNK")    
                info.SetResidueNumber(1)    
                info.SetChainId("X")    
                atom.SetPDBResidueInfo(info)    
        
        return rdmol    
    
    except Exception as e:    
        print(f"DEBUG: Exception during SDF processing: {e}")    
        return None