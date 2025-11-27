import os  
import numpy  
from IPython.core.debugger import set_trace  
from subprocess import Popen, PIPE  
import pymesh  
  
from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin  
import random  
  
# 导入新的混合方法  
try:  
    from triangulation.computeAPBS_openbabel import computeAPBS_hybrid, has_unsupported_atoms  
    OPENBABEL_AVAILABLE = True  
except ImportError:  
    print("警告: 无法导入OpenBabel混合方法,将仅使用标准PDB2PQR")  
    OPENBABEL_AVAILABLE = False  
  
"""  
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.  
Pablo Gainza - LPDI STI EPFL 2019  
This file is part of MaSIF.  
Released under an Apache License 2.0  
"""  
  
def computeAPBS(vertices, pdb_file, tmp_file_base, mol2_file=None):  
    """  
    Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex  
    如果检测到不支持的原子类型,自动切换到OpenBabel混合方法  
    """  
      
    # 检查是否需要使用混合方法  
    if OPENBABEL_AVAILABLE and mol2_file and has_unsupported_atoms(mol2_file):  
        print("检测到PDB2PQR不支持的原子类型,切换到OpenBabel混合方法")  
        try:  
            return computeAPBS_hybrid(vertices, pdb_file, tmp_file_base, mol2_file)  
        except Exception as e:  
            print(f"OpenBabel混合方法失败: {e}")  
            print("回退到标准PDB2PQR方法(可能会失败)")  
            # 继续执行标准流程  
      
    # 标准PDB2PQR流程  
    fields = tmp_file_base.split("/")[0:-1]  
    directory = "/".join(fields) + "/"  
    filename_base = tmp_file_base.split("/")[-1]  
    pdbname = pdb_file.split("/")[-1]  
      
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
  
    # 过滤多重占位警告,只显示真正的错误  
    stderr_text = stderr.decode('utf-8')  
    if 'Error' in stderr_text or 'Traceback' in stderr_text:  
        print("### PDB2PQR ###\n", stderr_text)  
    elif 'KeyError' in stderr_text:  
        print("### PDB2PQR 错误 ###")  
        print(stderr_text)  
        print("提示: 如果是不支持的原子类型,请确保已安装OpenBabel并创建computeAPBS_openbabel.py")  
  
    args = [apbs_bin, filename_base + ".in"]  
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
    stdout, stderr = p2.communicate()  
      
    vertfile = open(directory + "/" + filename_base + ".csv", "w")  
    for vert in vertices:  
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))  
    vertfile.close()  
  
    print("### APBS ###\n", stderr.decode('utf-8'))  
  
    args = [  
        multivalue_bin,  
        filename_base + ".csv",  
        filename_base + ".dx",  
        filename_base + "_out.csv",  
    ]  
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)  
    stdout, stderr = p2.communicate()  
  
    print("### MULTIVALUE ###\n", stderr.decode('utf-8'))  
  
    # Read the charge file  
    chargefile = open(tmp_file_base + "_out.csv")  
    charges = numpy.array([0.0] * len(vertices))  
    for ix, line in enumerate(chargefile.readlines()):  
        charges[ix] = float(line.split(",")[3])  
  
    # remove_fn = os.path.join(directory, filename_base)  
    # # os.remove(remove_fn)  
    # # os.remove(remove_fn+'.csv')  
    # # os.remove(remove_fn+'.dx')  
    # # os.remove(remove_fn+'.in')  
    # # os.remove(remove_fn+'_out.csv')  
  
    return charges