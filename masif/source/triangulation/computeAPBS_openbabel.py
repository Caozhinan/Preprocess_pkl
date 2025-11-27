import os  
import numpy as np  
from openbabel import openbabel as ob  
from subprocess import Popen, PIPE  
  
def has_unsupported_atoms(mol2_file):  
    """  
    检查MOL2文件中是否包含PDB2PQR不支持的原子类型  
    """  
    if mol2_file is None or not os.path.exists(mol2_file):  
        return False  
      
    # PDB2PQR PEOE不支持的原子类型  
    unsupported_atoms = ['B', 'V', 'Se', 'Si', 'As', 'Sb', 'Te']  
    unsupported_types = ['S.O']  # MOL2原子类型  
      
    try:  
        with open(mol2_file, 'r') as f:  
            content = f.read()  
            in_atom_section = False  
              
            for line in content.split('\n'):  
                if '@<TRIPOS>ATOM' in line:  
                    in_atom_section = True  
                    continue  
                elif '@<TRIPOS>' in line and in_atom_section:  
                    break  
                      
                if in_atom_section and line.strip():  
                    parts = line.split()  
                    if len(parts) >= 6:  
                        atom_type = parts[5]  # MOL2原子类型在第6列  
                        atom_name = parts[1]  # 原子名在第2列  
                          
                        # 检查原子类型  
                        if atom_type in unsupported_types:  
                            return True  
                        # 检查原子名的第一个字符  
                        for unsup in unsupported_atoms:  
                            if atom_name.startswith(unsup):  
                                return True  
        return False  
    except Exception as e:  
        print(f"检查MOL2文件时出错: {e}")  
        return False  
  
  
def compute_charges_with_openbabel(mol2_file, pdb_file, output_pqr):  
    """改进版:更准确地匹配配体原子"""  
    try:  
        # 1. 读取MOL2文件获取原子信息  
        obConversion = ob.OBConversion()  
        obConversion.SetInFormat("mol2")  
          
        ligand_mol = ob.OBMol()  
        if not obConversion.ReadFile(ligand_mol, mol2_file):  
            print(f"无法读取MOL2文件: {mol2_file}")  
            return False  
          
        # 计算Gasteiger电荷  
        charge_model = ob.OBChargeModel.FindType("gasteiger")  
        if not charge_model:  
            print("无法找到Gasteiger电荷模型")  
            return False  
        charge_model.ComputeCharges(ligand_mol)  
          
        # 2. 从MOL2文件读取原子名称和坐标  
        mol2_atoms = []  
        with open(mol2_file, 'r') as f:  
            in_atom_section = False  
            for line in f:  
                if '@<TRIPOS>ATOM' in line:  
                    in_atom_section = True  
                    continue  
                elif '@<TRIPOS>' in line and in_atom_section:  
                    break  
                if in_atom_section and line.strip():  
                    parts = line.split()  
                    if len(parts) >= 6:  
                        mol2_atoms.append({  
                            'name': parts[1],  
                            'x': float(parts[2]),  
                            'y': float(parts[3]),  
                            'z': float(parts[4]),  
                            'type': parts[5]  
                        })  
          
        # 3. 为每个原子生成PQR行  
        ligand_pqr_lines = []  
        for idx, atom_info in enumerate(mol2_atoms):  
            ob_atom = ligand_mol.GetAtom(idx + 1)  
            charge = ob_atom.GetPartialCharge()  
            radius = ob.GetVdwRad(ob_atom.GetAtomicNum())  
              
            # 使用MOL2中的坐标和原子名  
            pqr_line = f"HETATM{idx+1:5d}  {atom_info['name']:4s} LIG X   1    {atom_info['x']:8.3f}{atom_info['y']:8.3f}{atom_info['z']:8.3f} {charge:7.4f} {radius:6.4f}\n"  
            ligand_pqr_lines.append(pqr_line)  
          
        # 写入配体PQR文件  
        with open(output_pqr, 'w') as f:  
            f.writelines(ligand_pqr_lines)  
          
        print(f"成功生成配体PQR文件,包含{len(ligand_pqr_lines)}个原子")  
        return True  
          
    except Exception as e:  
        print(f"OpenBabel电荷计算失败: {e}")  
        import traceback  
        traceback.print_exc()  
        return False
  
def computeAPBS_hybrid(vertices, pdb_file, tmp_file_base, mol2_file=None):  
    """  
    混合方法:蛋白质用PDB2PQR,配体用OpenBabel  
    """  
    from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin  
      
    fields = tmp_file_base.split("/")[0:-1]  
    directory = "/".join(fields) + "/"  
    filename_base = tmp_file_base.split("/")[-1]  
    pdbname = pdb_file.split("/")[-1]  
      
    # 检查是否需要使用OpenBabel  
    use_openbabel = mol2_file and has_unsupported_atoms(mol2_file)  
      
    if use_openbabel:  
        print(f"检测到不支持的原子类型,使用OpenBabel计算配体电荷...")  
          
        # 1. 先用PDB2PQR处理蛋白质(不包含配体)  
        protein_only_pdb = os.path.join(directory, filename_base + "_protein.pdb")  
        with open(pdb_file, 'r') as fin, open(protein_only_pdb, 'w') as fout:  
            for line in fin:  
                if line.startswith('ATOM'):  # 只保留蛋白质原子  
                    fout.write(line)  
          
        # 运行PDB2PQR处理蛋白质  
        args = [  
            pdb2pqr_bin,  
            os.path.basename(protein_only_pdb),  
            filename_base + "_protein",  
            "--ff=PARSE",  
            "--whitespace",  
            "--noopt",  
            "--nodebump",  
            "--assign-only",  
            "--chain",  
            "--apbs-input",  
        ]  
          
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
        stdout, stderr = p2.communicate()  
          
        # 2. 用OpenBabel处理配体  
        ligand_pqr = os.path.join(directory, filename_base + "_ligand.pqr")  
        if not compute_charges_with_openbabel(mol2_file, pdb_file, ligand_pqr):  
            print("OpenBabel处理失败,回退到零电荷")  
            return np.zeros(len(vertices))  
          
        # 3. 合并蛋白质和配体的PQR文件  
        protein_pqr = os.path.join(directory, filename_base + "_protein.pqr")  
        combined_pqr = os.path.join(directory, filename_base + ".pqr")  
          
        with open(combined_pqr, 'w') as fout:  
            if os.path.exists(protein_pqr):  
                with open(protein_pqr, 'r') as fin:  
                    fout.write(fin.read())  
            if os.path.exists(ligand_pqr):  
                with open(ligand_pqr, 'r') as fin:  
                    fout.write(fin.read())  
          
        # 4. 手动创建APBS输入文件  
        apbs_input = os.path.join(directory, filename_base + ".in")  
        with open(apbs_input, 'w') as f:  
            f.write(f"""read  
    mol pqr {filename_base}.pqr  
end  
elec  
    mg-auto  
    dime 97 97 97  
    cglen 100 100 100  
    fglen 80 80 80  
    cgcent mol 1  
    fgcent mol 1  
    mol 1  
    lpbe  
    bcfl sdh  
    pdie 2.0  
    sdie 78.54  
    srfm smol  
    chgm spl2  
    sdens 10.0  
    srad 1.4  
    swin 0.3  
    temp 298.15  
    calcenergy total  
    calcforce no  
    write pot dx {filename_base}  
end  
quit  
""")  
          
    else:  
        # 使用标准PDB2PQR流程  
        args = [  
            pdb2pqr_bin, pdbname, filename_base,  
            "--ff=PARSE",  
            "--whitespace",  
            "--noopt",  
            "--nodebump",  
            "--assign-only",  
            "--chain",  
            "--apbs-input",  
        ]  
        if mol2_file is not None:  
            args.append("--ligand=" + mol2_file)  
          
        p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
        stdout, stderr = p2.communicate()  
      
    # 运行APBS  
    args = [apbs_bin, filename_base + ".in"]  
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
    stdout, stderr = p2.communicate()  
      
    # 写入顶点坐标  
    vertfile = open(directory + "/" + filename_base + ".csv", "w")  
    for vert in vertices:  
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))  
    vertfile.close()  
      
    print("### APBS ###\n", stderr.decode('utf-8'))  
      
    # 运行MULTIVALUE  
    args = [  
        multivalue_bin,  
        filename_base + ".csv",  
        filename_base + ".dx",  
        filename_base + "_out.csv",  
    ]  
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
    stdout, stderr = p2.communicate()  
      
    print("### MULTIVALUE ###\n", stderr.decode('utf-8'))  
      
    # 读取电荷文件  
    chargefile = open(tmp_file_base + "_out.csv")  
    charges = np.array([0.0] * len(vertices))  
    for ix, line in enumerate(chargefile.readlines()):  
        charges[ix] = float(line.split(",")[3])  
      
    # 清理临时文件  
    remove_fn = os.path.join(directory, filename_base)  
    for ext in ['', '.csv', '.dx', '.in', '_out.csv', '_protein.pdb', '_protein.pqr', '_ligand.pqr']:  
        try:  
            if os.path.exists(remove_fn + ext):  
                os.remove(remove_fn + ext)  
        except:  
            pass  
      
    return charges