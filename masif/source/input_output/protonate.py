"""
protonate.py: Wrapper method for the reduce program: protonate (i.e., add hydrogens) a pdb using reduce 
                and save to an output file.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""

from subprocess import Popen, PIPE
from IPython.core.debugger import set_trace
import os


def protonate(in_pdb_file, out_pdb_file):
    # protonate (i.e., add hydrogens) a pdb using reduce and save to an output file.
    # in_pdb_file: file to protonate.
    # out_pdb_file: output file where to save the protonated pdb file. 
    
    # Remove protons first, in case the structure is already protonated
    args = ["reduce", "-Trim", in_pdb_file]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()
    outfile = open(out_pdb_file, "w")
    outfile.write(stdout.decode('utf-8').rstrip())
    outfile.close()

    # Now add them again.
    # args = ["reduce", "-HIS", "-DB", "/work/upcorreia/bin/reduce/reduce_wwPDB_het_dict_old.txt", out_pdb_file]
    args = ["reduce", out_pdb_file, "-HIS"]
    het_dict = os.environ.get('REDUCE_HET_DICT')
    if het_dict is not None:
        args.extend(["-DB", het_dict])
    
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()
    outfile = open(out_pdb_file, "w")
    outfile.write(stdout.decode('utf-8'))
    outfile.close()

def fix_ligand_atom_names(pdb_file):  
    """修复配体中重复的原子名"""  
    with open(pdb_file, 'r') as f:  
        lines = f.readlines()  
      
    # 统计每个残基中的原子名出现次数  
    residue_atom_counts = {}  
    fixed_lines = []  
      
    for line in lines:  
        if line.startswith('HETATM') and 'UNK' in line:  
            # 解析PDB行  
            atom_name = line[12:16].strip()  
            res_name = line[17:20].strip()  
            chain_id = line[21:22].strip()  
            res_seq = line[22:26].strip()  
              
            residue_key = f"{res_name}_{chain_id}_{res_seq}"  
              
            if residue_key not in residue_atom_counts:  
                residue_atom_counts[residue_key] = {}  
              
            if atom_name in residue_atom_counts[residue_key]:  
                # 原子名重复，需要重命名  
                residue_atom_counts[residue_key][atom_name] += 1  
                new_atom_name = f"{atom_name}{residue_atom_counts[residue_key][atom_name]}"  
            else:  
                residue_atom_counts[residue_key][atom_name] = 1  
                new_atom_name = atom_name  
              
            # 重构PDB行，确保原子名4字符对齐  
            new_line = line[:12] + f"{new_atom_name:<4}" + line[16:]  
            fixed_lines.append(new_line)  
        else:  
            fixed_lines.append(line)  
      
    # 写回文件  
    with open(pdb_file, 'w') as f:  
        f.writelines(fixed_lines)  