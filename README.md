# __AffinSculptor 预处理脚本使用指南__

## 概述 
preprocess.sh 是一个自动化脚本，用于批量处理蛋白质-配体复合物数据，生成包含 MaSIF 表面指纹的整合特征文件。

## 前置要求 
环境配置：确保已安装所需的 Python 环境和依赖包
脚本权限：给脚本添加执行权限
chmod +x preprocess.sh
## 输入格式 
脚本需要一个 CSV 文件作为输入，格式如下：
receptor,ligand,name,pk,rmsd  
/xcfhome/zncao02/dataset_bap/test_set/custom/6uux-QHM/protein.pdb,/xcfhome/zncao02/dataset_bap/test_set/custom/6uux-QHM/ligand.sdf,6uux-QHM,6.63,1.25  
/xcfhome/zncao02/dataset_bap/test_set/custom/1abc-XYZ/protein.pdb,/xcfhome/zncao02/dataset_bap/test_set/custom/1abc-XYZ/ligand.sdf,1abc-XYZ,7.2,0.8
CSV 字段说明 
receptor: 蛋白质 PDB 文件的绝对路径
ligand: 配体 SDF 文件的绝对路径
name: 复合物名称标识符
pk: pK 值（数值）
rmsd: RMSD 值（数值）
使用方法 
基本调用 
./preprocess.sh input_data.csv
## 完整示例 
### 1. 准备输入 CSV 文件  
cat > my_complexes.csv << EOF  
receptor,ligand,name,pk,rmsd  
/path/to/complex1/protein.pdb,/path/to/complex1/ligand.sdf,complex1,6.5,1.2  
/path/to/complex2/protein.pdb,/path/to/complex2/ligand.sdf,complex2,7.1,0.9  
EOF  
  
### 2. 运行脚本  
./preprocess.sh my_complexes.csv
处理流程 
脚本会依次执行以下步骤：

Step 1: 调用 custom_input.py 进行预处理，生成基础 PKL 文件
Step 2: 调用 meshfeatureGen.py 生成网格和表面特征
Step 3: 调用 feature_precompute.py 进行特征预计算
Step 4: 调用 fingerprint_gen.py 生成 MaSIF 指纹
Step 5: 调用 merge_pkl.py 整合所有特征到最终 PKL 文件
输出结果 
对于每个复合物，脚本会在相应目录下创建 output 文件夹，包含：

/path/to/complex/output/  
├── complex.pdb                    # 预处理后的复合物结构  
├── surfaces/                      # 表面网格文件  
├── precomputed/                   # 预计算特征  
├── descriptors/                   # MaSIF 指纹  
├── original_data.pkl              # 原始图数据  
└── original_data_with_masif.pkl   # 整合 MaSIF 特征的最终文件  

## 错误处理 
如果某个复合物处理失败，脚本会继续处理下一个
每个步骤都有错误检查，失败时会显示具体错误信息
处理完成后会显示总体处理结果
## 注意事项 
路径要求：CSV 中的文件路径必须是绝对路径
文件存在性：确保所有输入的 PDB 和 SDF 文件都存在
磁盘空间：每个复合物会生成较多中间文件，确保有足够磁盘空间
处理时间：MaSIF 指纹生成可能需要较长时间，特别是对于大型复合物
## 故障排除 
常见错误 
权限错误：确保脚本有执行权限
路径错误：检查 CSV 中的文件路径是否正确
依赖缺失：确保所有 Python 脚本都已正确配置命令行参数支持
调试模式 
可以手动执行单个步骤来调试问题：

## 测试单个复合物的处理  
python /xcfhome/zncao02/AffinSculptor/preprocess/custom_input.py \\  
    --receptor /path/to/protein.pdb \\  
    --ligand /path/to/ligand.sdf \\  
    --output_dir /path/to/output \\  
    --name test_complex

这个脚本设计为批量处理工具，可以高效地为多个蛋白质-配体复合物生成包含 MaSIF 表面指纹的完整特征数据。