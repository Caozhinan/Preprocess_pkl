#!/bin/bash

source /es01/paratera/parasoft/module.sh
module load mambaforge/24.11.0-1-hxl

# 初始化 conda
eval "$(/es01/paratera/parasoft/soft/mambaforge/24.11.0-1/bin/conda shell.bash hook)"
conda activate /es01/paratera/sce0413/czn/conda_env/affincraft

# 设置环境变量
export TMPDIR="/tmp/masif_1p01_$$"
export MASIF_TMP_DIR="$TMPDIR"
mkdir -p "$TMPDIR"

# 运行 meshfeatureGen.py
python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/meshfeatureGen.py \
    --pdb_file "complex.pdb" \
    --chain_id "A" \
    --ligand_code "UNK" \
    --ligand_chain "X" \
    --sdf_file "6rsa_ligand.sdf" \
    --output_dir "output"

# 清理临时目录
rm -rf "$TMPDIR"