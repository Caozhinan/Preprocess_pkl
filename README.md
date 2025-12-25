# Preprocess_pkl: Protein-Ligand Complex Feature Extraction Pipeline 
A high-performance parallel preprocessing pipeline for extracting and integrating structural features from protein-ligand complexes, combining traditional molecular features with MaSIF surface fingerprints for machine learning applications.

## Overview 
This repository implements a comprehensive preprocessing pipeline that transforms protein-ligand complexes into feature-rich PKL files suitable for deep learning models. The pipeline integrates molecular graph features with MaSIF (Machine Learning for Structural Interaction Fingerprinting) surface descriptors through a multi-step parallel processing workflow .

## Key Features 
Parallel Processing: Utilizes GNU Parallel for efficient multi-core processing of large datasets
Integrated Features: Combines traditional molecular features with MaSIF surface fingerprints
Robust Error Handling: Continues processing on individual failures with comprehensive logging
Resource Isolation: Each worker uses isolated temporary directories to prevent conflicts
Flexible Input: Supports both PDB and SDF ligand formats custom_input.py:35-56
Pipeline Architecture 
The processing pipeline consists of four main steps:







### Step 1: Base Feature Extraction 
Generates molecular graph features using custom_input.py, including coordinates, node features, edge indices, and chemical properties .

### Step 2: Surface Mesh Generation 
Creates surface meshes for both ligand and pocket using MSMS algorithm with configurable density and resolution parameters .

### Step 3: MaSIF Feature Computation 
Computes geometric and chemical surface features using MaSIF's precomputation pipeline with configurable distance and shape parameters .

### Step 4: Feature Integration 
Merges base molecular features with MaSIF surface descriptors into unified PKL files .

Installation 
Prerequisites 
Linux environment with GNU Parallel
Conda environment with required packages
Python 3.7+ with scientific computing libraries
Environment Setup 
# Activate conda environment (adjust path as needed)  
conda activate /path/to/your/conda_env/affincraft
Usage 
Basic Usage 
./update_preprocess.sh <csv_file> [num_cores] [log_dir]
Parameters 
csv_file: Path to CSV file with format receptor,ligand,name,pk,rmsd
num_cores: Number of parallel cores (default: 1)
log_dir: Custom log directory (default: ./logs/<csv_name>)
Example 
## Process dataset using 8 cores  
./update_preprocess.sh complexes.csv 8 ./processing_logs  
  
## Single-core processing with default logging  
./update_preprocess.sh complexes.csv
Input Format 
The CSV file should contain the following columns:

receptor,ligand,name,pk,rmsd  
/path/to/receptor.pdb,/path/to/ligand.sdf,complex_name,1.5,0.8  
Output Structure 
For each complex, the pipeline generates:

{name}_features.pkl: Base molecular features
{name}_features_with_masif.pkl: Integrated features with MaSIF descriptors
Intermediate files are automatically cleaned up after successful processing .

Logging and Monitoring 
The pipeline provides comprehensive logging:

<<<<<<< HEAD
## 测试单个复合物的处理  
python /xcfhome/zncao02/AffinSculptor/preprocess/custom_input.py \\  
    --receptor /path/to/protein.pdb \\  
    --ligand /path/to/ligand.sdf \\  
    --output_dir /path/to/output \\  
    --name test_complex

这个脚本设计为批量处理工具，可以高效地为多个蛋白质-配体复合物生成包含 MaSIF 表面指纹的完整特征数据。


## 关键技术更新 
1. 输入处理优化
- 修改 custom_input.py: 根据配体文件后缀(.sdf或.pdb)自动选择处理函数
- 配体残基名标准化: 将配体残基名从UNK改为LIG preprocess.py:58-67
- 溶剂过滤: 自动过滤常见溶剂分子，避免PLIP误识别 preprocess.py:58-67
2. 小分子表面生成技术突破
- MSMS参数优化: density=8.0, probe_radius=1.2Å, mesh_res=0.8
- 静电特征解决方案: 使用Gasteiger电荷替代APBS，解决小分子APBS计算失败问题
- XYZRN转换修复: 解决小分子原子过滤问题，确保MSMS正常运行
3. 特征预计算参数调整
- 小分子参数: max_distance=5.0Å, max_shape_size=60
- 蛋白质口袋参数: max_distance=6.0Å, max_shape_size=80
##  新增脚本列表 
meshfeatureGen_ligand.py	小分子表面生成	高密度MSMS参数，Gasteiger电荷
meshfeatureGen_pocket.py	蛋白质口袋表面生成	跳过链ID检测，直接处理所有残基
feature_precompute_ligand.py	小分子特征预计算	5维特征，184个patches
feature_precompute_pocket.py	蛋白质口袋特征预计算	5维特征，1759个patches
merge_surface_features.py	特征合并	过滤冗余字段，添加MaSIF特征
batch_preprocess_surface_features_full.sh	完整批量处理	四步流程，多核并行，自动清理
## 输出文件结构 
每个复合物最终生成：
- output/precomputed/ligand/all_features.npz - 小分子5维特征 (184×60×5)
- output/precomputed/pocket/all_features.npz - 蛋白质口袋5维特征 (1759×80×5)
- output/{name}_features_with_masif.pkl - 合并后的完整特征文件
=======
parallel.log: GNU Parallel job tracking
output_success.log: Success messages and progress
output_fail.log: Error messages and failures
Performance Considerations 
Memory: Each worker typically requires 2-4 GB RAM
Disk: Temporary files during processing (~500-1000 MB per complex)
CPU: Steps 3-4 are CPU-intensive; Steps 1-2 are I/O-bound
Recommended Cores: Use min(CPU_CORES - 2, MEMORY_GB / 4) for optimal performance
>>>>>>> 4f9f0f7 (ppi_to_be_done)
