#!/bin/bash



if [ $# -lt 1 ]; then
    echo "Usage: $0 <csv_file> [merged_pkl_path] [num_cores]"
    echo "CSV format: receptor,ligand,name,pk,rmsd"
    exit 1
fi

CSV_FILE="$1"
MERGED_PKL="$2"
NUM_CORES="${3:-1}"

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

    export TMPDIR="$PROCESS_TMP_DIR"
    export MASIF_TMP_DIR="$PROCESS_TMP_DIR"

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
    python3 /xcfhome/zncao02/affinsculp/preprocess/custom_input.py \
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
    python3 /xcfhome/zncao02/affinsculp/masif/meshfeatureGen.py \
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
    python3 /xcfhome/zncao02/affinsculp/masif/feature_precompute.py \
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
    python3 /xcfhome/zncao02/affinsculp/masif/fingerprint_gen.py \
        --precomputed_dir "$precomputed_dir" \
        --output_dir "$output_dir" \
        --ppi_pair_id "$ppi_pair_id" \
        --custom_params_file "/xcfhome/zncao02/affinsculp/masif/data/masif_ppi_search/nn_models/sc05/all_feat/custom_params.py" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in fingerprint_gen.py for $name"
        return 1
    fi

    echo "Step 5: Merging MaSIF features to existing PKL..."
    python3 /xcfhome/zncao02/affinsculp/preprocess/merge_pkl.py \
        --pkl_file "$pkl_output" \
        --output_dir "$output_dir" \
        --name "$name" \
        --pk "$pk" \
        --rmsd "$rmsd" \
        --chain_id "$protein_chain" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in merge_pkl.py for $name"
        return 1
    fi

    echo "Successfully processed: $name"
    return 0
}

export -f process_one_complex
export -f detect_protein_chain

# total_lines=$(tail -n +2 "$CSV_FILE" | wc -l)
# echo "Total complexes to process: $total_lines"

# 创建日志文件夹
LOG_DIR="/xcfhome/zncao02/affinsculp/log_temp/logs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOG_DIR"

# 正确分隔符传递到parallel（用tab）
tail -n +2 "$CSV_FILE" | \
    awk -F',' '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | \
    parallel --colsep '\t' -j "$NUM_CORES" --joblog "$LOG_DIR/parallel.log" \
    process_one_complex {1} {2} {3} {4} {5} \
    >"$LOG_DIR/output_success.log" 2>"$LOG_DIR/output_fail.log"

# 统计成功/失败
success_count=$(grep -c "Successfully processed:" "$LOG_DIR/output_success.log" 2>/dev/null || echo "0")
fail_count=$(grep -c "Error" "$LOG_DIR/output_fail.log" 2>/dev/null || echo "0")

echo "Processing finished! Success: $success_count, Failed: $fail_count"
echo "Logs saved in: $LOG_DIR"

echo "All processing completed!"