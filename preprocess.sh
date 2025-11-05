#!/bin/bash
# source /es01/paratera/parasoft/module.sh
# module load mambaforge/24.11.0-1-hxl
source /es01/paratera/parasoft/soft/mambaforge/24.11.0-1/etc/profile.d/conda.sh
# 初始化 conda
# eval "$(conda shell.bash hook)"
# 激活环境
conda activate /es01/paratera/sce0413/czn/conda_env/affincraft
# which python

if [ $# -lt 1 ]; then
    echo "Usage: $0 <csv_file> [num_cores]"
    echo "CSV format: receptor,ligand,name,pk,rmsd"
    exit 1
fi

CSV_FILE="$1"
NUM_CORES="${2:-1}"

if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file $CSV_FILE not found"
    exit 1
fi

detect_protein_chain() {
    local pdb_file="$1"
    local protein_chain
    protein_chain=$(grep "^ATOM" "$pdb_file" | awk '{print $5}' | grep -v "X" | sort | uniq -c | sort -nr | head -1 | awk '{print $2}')
    echo "$protein_chain"
}

process_one_complex() {
    local receptor="$1"
    local ligand="$2"
    local name="$3"
    local pk="$4"
    local rmsd="$5"

    [ -z "$receptor" ] && return 0

    echo "Processing: $name"
    target_dir=$(dirname "$receptor")
    output_dir="$target_dir/output"
    mkdir -p "$output_dir"

    THREAD_ID="${name}_$(date +%s%N)_$$"
    PROCESS_TMP_DIR="/tmp/masif_${THREAD_ID}"
    mkdir -p "$PROCESS_TMP_DIR"

    cleanup() {
        rm -rf "$PROCESS_TMP_DIR" 2>/dev/null || true
        rm -f "$temp_csv" 2>/dev/null || true
    }
    trap cleanup EXIT

    temp_csv="$output_dir/temp_input_${THREAD_ID}.csv"
    echo "receptor,ligand,name,pk,rmsd" > "$temp_csv"
    echo "$receptor,$ligand,$name,$pk,$rmsd" >> "$temp_csv"

    pkl_output="$output_dir/${name}_features.pkl"

    echo "Step 1: Running custom_input.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/preprocess/custom_input.py \
        "$temp_csv" \
        "$pkl_output" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in custom_input.py for $name"
        return 1
    fi
    rm -f "$temp_csv"

    if [ ! -f "$pkl_output" ]; then
        echo "Error: PKL file not generated: $pkl_output"
        return 1
    fi

    complex_pdb="$target_dir/complex.pdb"
    if [ ! -f "$complex_pdb" ]; then
        echo "Error: Complex PDB file not found: $complex_pdb"
        return 1
    fi

    protein_chain=$(detect_protein_chain "$complex_pdb")
    if [ -z "$protein_chain" ]; then
        echo "Error: Could not detect protein chain in $complex_pdb"
        return 1
    fi

    echo "Step 2: Running meshfeatureGen.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/meshfeatureGen.py \
        --pdb_file "$complex_pdb" \
        --chain_id "$protein_chain" \
        --ligand_code "UNK" \
        --ligand_chain "X" \
        --sdf_file "$ligand" \
        --output_dir "$output_dir" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in meshfeatureGen.py for $name"
        return 1
    fi

    ply_file="$output_dir/surfaces/complex_${protein_chain}.ply"
    precomputed_dir="$output_dir/precomputed/complex_${protein_chain}"
    ppi_pair_id="complex_${protein_chain}"

    if [ ! -f "$ply_file" ]; then
        echo "Error: PLY file not generated: $ply_file"
        return 1
    fi

    echo "Step 3: Running feature_precompute.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/feature_precompute.py \
        --ply_file "$ply_file" \
        --output_dir "$output_dir" \
        --masif_app "masif_site" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in feature_precompute.py for $name"
        return 1
    fi

    if [ ! -d "$precomputed_dir" ]; then
        echo "Error: Precomputed directory not generated: $precomputed_dir"
        return 1
    fi

    echo "Step 4: Running fingerprint_gen.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/fingerprint_gen.py \
        --precomputed_dir "$precomputed_dir" \
        --output_dir "$output_dir" \
        --ppi_pair_id "$ppi_pair_id" \
        --custom_params_file "/es01/paratera/sce0413/czn/preprocess_pkl/masif/data/masif_ppi_search/nn_models/sc05/all_feat/custom_params.py" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in fingerprint_gen.py for $name"
        return 1
    fi

    echo "Step 5: Merging MaSIF features to existing PKL..."
    descriptors_dir="$output_dir/descriptors/complex_${protein_chain}"

    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/preprocess/merge_pkl.py \
        --pkl_file "$pkl_output" \
        --descriptors_dir "$descriptors_dir" \
        --output_file "$output_dir/${name}_features_with_masif.pkl" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in merge_pkl.py for $name"
        return 1
    fi

    echo "Successfully processed: $name"
    return 0
}

export -f process_one_complex
export -f detect_protein_chain

# 创建日志文件夹
LOG_DIR="/es01/paratera/sce0413/czn/preprocess_pkl/log-part4/$(basename "$CSV_FILE" .csv)"
mkdir -p "$LOG_DIR"

# 正确分隔符传递到parallel（tab分隔，并去掉空格）
tail -n +2 "$CSV_FILE" | \
    awk -F',' '{
        gsub(/^ +| +$/, "", $1);
        gsub(/^ +| +$/, "", $2);
        gsub(/^ +| +$/, "", $3);
        gsub(/^ +| +$/, "", $4);
        gsub(/^ +| +$/, "", $5);
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5
    }' | \
    parallel --colsep '\t' -j "$NUM_CORES" --joblog "$LOG_DIR/parallel.log" \
    process_one_complex {1} {2} {3} {4} {5} \
    >"$LOG_DIR/output_success.log" 2>"$LOG_DIR/output_fail.log"

# 统计成功/失败
success_count=$(grep -c "Successfully processed:" "$LOG_DIR/output_success.log" 2>/dev/null || echo "0")
fail_count=$(grep -c "Error" "$LOG_DIR/output_fail.log" 2>/dev/null || echo "0")

echo "Processing finished! Success: $success_count, Failed: $fail_count"
echo "Logs saved in: $LOG_DIR"

if [ $success_count -gt 0 ]; then
    echo "Starting cleanup and merge process..."

    # 创建合并输出目录
    MERGE_OUTPUT_DIR="/es01/paratera/sce0413/czn/preprocess_pkl/merged_results4/$(basename "$CSV_FILE" .csv)"
    mkdir -p "$MERGE_OUTPUT_DIR"

    # 运行清理和合并脚本（注意：反斜杠后不能有空格）
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/preprocess/merge_all_pkl.py \
        --csv_file "$CSV_FILE" \
        --output_dir "$MERGE_OUTPUT_DIR"

    if [ $? -eq 0 ]; then
        echo "Cleanup and merge completed successfully!"
        echo "Merged files saved in: $MERGE_OUTPUT_DIR"
    else
        echo "Error in cleanup and merge process"
    fi
else
    echo "No successful processing found, skipping cleanup and merge"
fi

echo "All processing completed!"