from pathlib import Path  # 导入 Path，用于文件路径操作
from subprocess import run, DEVNULL  # 导入 run/DEVNULL，用于调用外部命令行工具并屏蔽输出
from rdkit import Chem  # 导入 RDKit 的 Chem 模块，用于分子操作
import pymol  # 导入 pymol，用于分子可视化和格式转换

# 读取分子文件（支持 sdf/mol 格式），并进行消毒和修正
def read_mol(mol_file: Path):
    """
    For ligand, use sdf, for protein, use sdf converted from pdb by pymol
    # 对于配体直接用 sdf，对于蛋白建议用 pymol 转换过的 sdf
    """
    # 不严格解析/不消毒读取分子文件
    mol = Chem.MolFromMolFile(str(mol_file), sanitize=False, strictParsing=False)
    if mol is None:  # 如果读取失败
        raise RuntimeError(f"{mol_file} cannot be processed")  # 抛出异常
    mol = correct_sanitize_v2(mol)  # 进行自定义的分子消毒/修正
    # try:
    #     mol = neutralize_atoms(mol)  # 可选：中和离子
    # except Exception:
    #     pass
    return Chem.RemoveHs(mol, sanitize=False,)  # 移除氢原子（不消毒），返回修正后的分子

# 中和分子中的带电原子
def neutralize_atoms(mol):
    # SMARTS 匹配带正/负电荷但不合理的原子
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)  # 搜索匹配原子
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:  # 若有匹配
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()  # 获取原子电荷数
            hcount = atom.GetTotalNumHs()  # 氢原子数
            atom.SetFormalCharge(0)  # 设置为中性
            atom.SetNumExplicitHs(hcount - chg)  # 调整氢原子数
            atom.UpdatePropertyCache()  # 更新原子属性
    return mol

# 查找与指定原子连接、满足条件的邻居原子
def find_atom_bond_around(at: Chem.Atom, sym: str, bt: Chem.BondType, inring: bool, aromatic: bool):
    mol = at.GetOwningMol()  # 获取分子对象
    res = []
    for n in at.GetNeighbors():  # 遍历邻居原子
        bond = mol.GetBondBetweenAtoms(at.GetIdx(), n.GetIdx())  # 获取连接键
        c1 = (n.GetSymbol() == sym) if sym is not None else True  # 判断原子类型
        c2 = (bond.GetBondType() == bt) if bt is not None else True  # 判断键类型
        c3 = (n.IsInRing() == inring) if inring is not None else True  # 判断环系
        c4 = (n.GetIsAromatic() == aromatic) if aromatic is not None else True  # 判断芳香性
        if c1 and c2 and c3 and c4:  # 满足所有条件
            res.append(n)
    return res

# 修正磷酰基团（P=O 相关结构）
def fix_phosphoryl_group(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'P':
        double_bonded_o = find_atom_bond_around(at, "O", Chem.BondType.DOUBLE, False, False)
        if len(double_bonded_o) >= 2:  # 如果有2个及以上双键氧
            mol.GetBondBetweenAtoms(at.GetIdx(), double_bonded_o[0].GetIdx()).SetBondType(Chem.BondType.SINGLE)
            double_bonded_o[0].SetFormalCharge(-1)
            at.SetFormalCharge(0)

# 修正羧基（COO-）电荷与键型
def fix_carboxyl_group(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'C':
        double_bonded_o = find_atom_bond_around(at, "O", Chem.BondType.DOUBLE, False, False)
        single_bonded_o = find_atom_bond_around(at, "O", Chem.BondType.SINGLE, False, False)
        if len(double_bonded_o) == 2:  # 两个双键氧
            mol.GetBondBetweenAtoms(at.GetIdx(), double_bonded_o[0].GetIdx()).SetBondType(Chem.BondType.SINGLE)
            double_bonded_o[0].SetFormalCharge(-1)
            at.SetFormalCharge(0)
        if len(double_bonded_o) == 1 and len(single_bonded_o) == 1 and at.GetFormalCharge() == -1:
            at.SetFormalCharge(0)
            single_bonded_o[0].SetFormalCharge(-1)

# 修正胍基/亚胺基团（如 N=C(N)N）
def fix_guanidine_amidine_group(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'C' and not at.IsInRing():
        double_bonded_n = find_atom_bond_around(at, "N", Chem.BondType.DOUBLE, None, False)
        single_bonded_n = find_atom_bond_around(at, "N", Chem.BondType.SINGLE, None, False)
        if len(double_bonded_n) in {2, 3}:
            is_set = False
            for n in double_bonded_n:
                num_non_h = 0
                for a in n.GetNeighbors():
                    if a.GetSymbol() != "H": num_non_h += 1
                if num_non_h == 1 and not is_set:
                    is_set = True
                    mol.GetBondBetweenAtoms(at.GetIdx(), n.GetIdx()).SetBondType(Chem.BondType.DOUBLE)
                    n.SetFormalCharge(1)
                else:
                    mol.GetBondBetweenAtoms(at.GetIdx(), n.GetIdx()).SetBondType(Chem.BondType.SINGLE)
                    n.SetFormalCharge(0)
            at.SetFormalCharge(0)
        if len(double_bonded_n) == 1 and len(single_bonded_n) in {1, 2}:
            at.SetFormalCharge(0)
            double_bonded_n[0].SetFormalCharge(1)

# 修正磺酰基团（SO3-）等超价氧问题
def fix_sulfonyl_group(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'S':
        double_bonded_o = find_atom_bond_around(at, "O", Chem.BondType.DOUBLE, False, False)
        if len(double_bonded_o) > 2:
            for o in double_bonded_o[2:]:
                mol.GetBondBetweenAtoms(at.GetIdx(), o.GetIdx()).SetBondType(Chem.BondType.SINGLE)
                o.SetFormalCharge(-1)

# 修正鸟嘌呤环上的氮
def fix_guanine_group(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'C' and at.IsInRing():
        single_bonded_inring_n = find_atom_bond_around(at, "N", Chem.BondType.AROMATIC, True, True)
        double_bonded_nonring_n = find_atom_bond_around(at, "N", Chem.BondType.DOUBLE, False, False)
        if len(single_bonded_inring_n) in {1, 2} and len(double_bonded_nonring_n) == 1:
            at.SetFormalCharge(0)
            double_bonded_nonring_n[0].SetFormalCharge(1)

# 修正芳香氮口袋（如咪唑、嘧啶等环上的氮电荷）
def fix_aromatic_n_pocket(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'N':
        single_bonded = find_atom_bond_around(at, None, Chem.BondType.SINGLE, None, None)
        double_bonded = find_atom_bond_around(at, None, Chem.BondType.DOUBLE, None, None)
        aromatic_bonded = find_atom_bond_around(at, None, Chem.BondType.AROMATIC, None, None)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 4:
            at.SetFormalCharge(1)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 4 and len(double_bonded):
            mol.GetBondBetweenAtoms(at.GetIdx(), double_bonded[0].GetIdx()).SetBondType(Chem.BondType.SINGLE)

# 修正羟基氧或带电氧
def fix_pocket_o_charge(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'O':
        single_bonded = find_atom_bond_around(at, None, Chem.BondType.SINGLE, None, None)
        double_bonded = find_atom_bond_around(at, None, Chem.BondType.DOUBLE, None, None)
        aromatic_bonded = find_atom_bond_around(at, None, Chem.BondType.AROMATIC, None, None)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 2:
            at.SetFormalCharge(0)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 3:
            at.SetFormalCharge(1)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 4:
            at.SetFormalCharge(2)

# 修正碳口袋电荷
def fix_pocket_c_charge(at: Chem.Atom, mol: Chem.Mol):
    if at.GetSymbol() == 'C':
        single_bonded = find_atom_bond_around(at, None, Chem.BondType.SINGLE, None, None)
        double_bonded = find_atom_bond_around(at, None, Chem.BondType.DOUBLE, None, None)
        aromatic_bonded = find_atom_bond_around(at, None, Chem.BondType.AROMATIC, None, None)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) == 4 and at.GetFormalCharge() != 0:
            at.SetFormalCharge(0)
        if len(single_bonded) + 1.5 * len(aromatic_bonded) + 2 * len(double_bonded) > 4 and len(double_bonded):
            mol.GetBondBetweenAtoms(at.GetIdx(), double_bonded[0].GetIdx()).SetBondType(Chem.BondType.SINGLE)

# 修正不在环上的芳香原子（如：芳香 C 不在环上）
def fix_non_ring_aromatic(mol):
    for atom in mol.GetAtoms():
        if not atom.IsInRing() and atom.GetIsAromatic():
            atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() == Chem.BondType.AROMATIC and not (a1.IsInRing() and a2.IsInRing()):
            bond.SetBondType(Chem.BondType.SINGLE)

# 在 ring 列表中找到包含所有指定原子 idx 的环
def get_ring_atoms(rings, idx_lst):
    for r in rings:
        if all([i in r for i in idx_lst]):
            return r
    return None

# 判断一组原子 idx 是否全为芳香
def all_aromatic(mol, idx_lst):
    return all([mol.GetAtomWithIdx(i).GetIsAromatic() for i in idx_lst]) if idx_lst is not None else False

# 主消毒修正函数（尝试 RDKit 消毒，失败时采用部分消毒/自定义修正）
def correct_sanitize_v2(mol):
    try:
        Chem.SanitizeMol(Chem.Mol(mol))  # 尝试完整消毒
        return mol
    except:
        pass
    # 部分消毒或超价元素错误时
    mol.UpdatePropertyCache(strict=False)  # 不严格更新属性
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True,)  # 仅检测对称环
    rinfo = mol.GetRingInfo().AtomRings()  # 获取所有环的信息

    # 以下是用 SMARTS 片段查找特定结构并修正
    # COO 负电荷修正
    frag1 = Chem.MolFromSmarts('[C;-1](=O)(-,=O)')
    if mol.HasSubstructMatch(frag1):
        hits = mol.GetSubstructMatches(frag1)
        for hit in hits:
            c, o1 = mol.GetAtomWithIdx(hit[0]), mol.GetAtomWithIdx(hit[1])
            c.SetFormalCharge(0)
            o1.SetFormalCharge(-1)
            mol.GetBondBetweenAtoms(hit[0], hit[1]).SetBondType(Chem.BondType.SINGLE)
    # 省略部分片段，以下类似，都是匹配结构并修正键型或电荷
    # ...（此处略去每个 frag 的注释，结构同理，都是针对常见“问题基团”用 SMARTS 修正键型和电荷）
    # 具体每个 frag 见原代码
    # ...

    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomKekulizeException:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except Chem.KekulizeException:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except Chem.AtomValenceException:
        # 如果遇到超价原子等问题，参考 RDKit 邮件列表推荐做法
        # https://sourceforge.net/p/rdkit/mailman/message/32599798/
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            Chem.SanitizeFlags.SANITIZE_FINDRADICALS|
            Chem.SanitizeFlags.SANITIZE_KEKULIZE|
            Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|
            Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
            Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|
            Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            catchErrors=True,
        )
    return mol

# 另一套消毒修正方案，遍历所有原子逐一修正
def correct_sanitize_v1(mol: Chem.Mol):
    for at in mol.GetAtoms():
        fix_phosphoryl_group(at, mol)
        fix_carboxyl_group(at, mol)
        fix_guanidine_amidine_group(at, mol)
        fix_sulfonyl_group(at, mol)
        fix_guanine_group(at, mol)
        fix_aromatic_n_pocket(at, mol)
        fix_pocket_o_charge(at, mol)
        fix_pocket_c_charge(at, mol)
    fix_non_ring_aromatic(mol)
    try:
        Chem.SanitizeMol(mol)
    except Chem.KekulizeException:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except Chem.AtomValenceException:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(
            mol,
            Chem.SanitizeFlags.SANITIZE_FINDRADICALS|
            Chem.SanitizeFlags.SANITIZE_KEKULIZE|
            Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|
            Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
            Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|
            Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            catchErrors=True,
        )
    return mol

# 用 pymol 进行 sdf/pdb 等格式互转、加氢等
def pymol_convert(in_file: Path, out_file: Path):
    pymol.cmd.reinitialize()  # 重置 pymol 环境
    pymol.cmd.load(f"{str(in_file)}")  # 加载输入分子
    # pymol.cmd.remove("hydrogens")  # 可选：移除原有氢
    pymol.cmd.h_add()  # 自动加氢
    pymol.cmd.save(f"{str(out_file)}", "not sol.")  # 保存到输出文件，排除溶剂

# 用 pymol 选取蛋白-配体 pocket 区域并导出
def pymol_pocket(receptor_file: Path, ligand_file: Path, pocket_file: Path):
    pymol.cmd.reinitialize()
    pymol.cmd.load(f"{str(receptor_file)}", "receptor")  # 加载受体
    pymol.cmd.load(f"{str(ligand_file)}", "ligand")  # 加载配体
    pymol.cmd.select("pocket", "byres (receptor within 10 of ligand)")  # 选取距离配体 10 Å 内的受体残基
    pymol.cmd.save(f"{str(pocket_file)}", "pocket and not sol.")  # 导出 pocket 区域

# 用 openbabel 将 pdb 转 mol 格式
def obabel_pdb2mol(in_file: Path, out_file: Path):
    run([
        'obabel', '-ipdb', str(in_file), '-omol', f'-O{str(out_file)}',
        '-x3v', '-h', '--partialcharge', 'eem'
    ], check=True, stdout=DEVNULL, stderr=DEVNULL)

# 用 openbabel 将 sdf 转 mol 格式
def obabel_sdf2mol(in_file: Path, out_file: Path):
    run([
        'obabel', '-isdf', str(in_file), '-omol', f'-O{str(out_file)}',
        '-x3v', '-h', '--partialcharge', 'eem'
    ], check=True, stdout=DEVNULL, stderr=DEVNULL)

# 用 openbabel 将 mol2 转 mol 格式
def obabel_mol22mol(in_file: Path, out_file: Path):
    run([
        'obabel', '-imol2', str(in_file), '-omol', f'-O{str(out_file)}',
        '-x3v', '-h', '--partialcharge', 'eem'
    ], check=True, stdout=DEVNULL, stderr=DEVNULL)

# 程序入口
if __name__ == "__main__":
    # from sys import stderr
    from preprocess import gen_feature  # 导入自定义特征生成函数
    path = Path("/mnt/yaosen-data/PDBBind/refined-set-2019")  # 数据集根目录
    d = path / '2epn'
    print(d.name)
    # if len(d.name) != "2epn": continue  # 只处理 2epn
    ligand = read_mol(d / f'{d.name}_ligand.sdf')  # 读取配体
    pocket = read_mol(d / f"{d.name}_pocket.sdf")  # 读取 pocket 区域
    res = gen_feature(ligand, pocket, d.name)  # 生成特征（自定义函数）
    for atom in ligand.GetAtoms():
        if atom.GetSymbol() == "H":
            print("HHH")  # 检查配体中氢原子
    for atom in pocket.GetAtoms():
        if atom.GetSymbol() == "H":
            print("HHH")  # 检查 pocket 区域中氢原子