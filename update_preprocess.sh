#!/bin/bash
# source /es01/paratera/parasoft/soft/mambaforge/24.11.0-1/etc/profile.d/conda.sh
# conda activate /es01/paratera/sce0413/czn/conda_env/affincraft

if [ $# -lt 1 ]; then
<<<<<<< HEAD
    echo "Usage: $0 <csv_file> [num_cores]"
=======
    echo "Usage: $0 <csv_file> [num_cores] [log_dir]"
>>>>>>> 4f9f0f7 (ppi_to_be_done)
    echo "CSV format: receptor,ligand,name,pk,rmsd"
    exit 1
fi

CSV_FILE="$1"
NUM_CORES="${2:-1}"
<<<<<<< HEAD
=======
LOG_DIR_ARG="${3:-}"
>>>>>>> 4f9f0f7 (ppi_to_be_done)

if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file $CSV_FILE not found"
    exit 1
fi

<<<<<<< HEAD
=======
# -------- log dir: now configurable by argument --------
if [ -n "$LOG_DIR_ARG" ]; then
    LOG_DIR="$LOG_DIR_ARG"
else
    # default if not provided (you can change this default)
    LOG_DIR="./logs/$(basename "$CSV_FILE" .csv)"
fi
mkdir -p "$LOG_DIR"

>>>>>>> 4f9f0f7 (ppi_to_be_done)
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

    # ========== Step 1: custom_input.py ==========
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

    pocket_pdb="$target_dir/${name}_protein_pocket.pdb"
    if [ ! -f "$pocket_pdb" ]; then
        echo "Error: Pocket PDB file not found: $pocket_pdb"
        return 1
    fi

    # ========== Step 2a: meshfeatureGen_ligand.py ==========
    echo "Step 2a: Running meshfeatureGen_ligand.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/meshfeatureGen_ligand.py \
        --sdf_file "$ligand" \
        --output_dir "$output_dir" \
        --density 8.0 \
        --probe 1.2 \
        --mesh_res 0.8 < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in meshfeatureGen_ligand.py for $name"
        return 1
    fi

    ligand_ply="$output_dir/surfaces/ligand.ply"
    if [ ! -f "$ligand_ply" ]; then
        echo "Error: Ligand PLY file not generated: $ligand_ply"
        return 1
    fi

    # ========== Step 2b: meshfeatureGen_pocket.py ==========
    echo "Step 2b: Running meshfeatureGen_pocket.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/meshfeatureGen_pocket.py \
        --pocket_pdb "$pocket_pdb" \
        --output_dir "$output_dir" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in meshfeatureGen_pocket.py for $name"
        return 1
    fi

    pocket_ply="$output_dir/surfaces/pocket.ply"
    if [ ! -f "$pocket_ply" ]; then
        echo "Error: Pocket PLY file not generated: $pocket_ply"
        return 1
    fi

    # ========== Step 3a: feature_precompute_ligand.py ==========
    echo "Step 3a: Running feature_precompute_ligand.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/feature_precompute_ligand.py \
        --ply_file "$ligand_ply" \
        --output_dir "$output_dir" \
        --max_distance 5.0 \
<<<<<<< HEAD
        --max_shape_size 60 < /dev/null
=======
        --max_shape_size 40 < /dev/null
>>>>>>> 4f9f0f7 (ppi_to_be_done)
    if [ $? -ne 0 ]; then
        echo "Error in feature_precompute_ligand.py for $name"
        return 1
    fi

    ligand_npz="$output_dir/precomputed/ligand/all_features.npz"
    if [ ! -f "$ligand_npz" ]; then
        echo "Error: Ligand features not generated: $ligand_npz"
        return 1
    fi

    # ========== Step 3b: feature_precompute_pocket.py ==========
    echo "Step 3b: Running feature_precompute_pocket.py..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/masif/feature_precompute_pocket.py \
        --ply_file "$pocket_ply" \
        --output_dir "$output_dir" \
        --max_distance 6.0 \
<<<<<<< HEAD
        --max_shape_size 60 < /dev/null
=======
        --max_shape_size 40 < /dev/null
>>>>>>> 4f9f0f7 (ppi_to_be_done)
    if [ $? -ne 0 ]; then
        echo "Error in feature_precompute_pocket.py for $name"
        return 1
    fi

    pocket_npz="$output_dir/precomputed/pocket/all_features.npz"
    if [ ! -f "$pocket_npz" ]; then
        echo "Error: Pocket features not generated: $pocket_npz"
        return 1
    fi

    # ========== Step 4: merge_surface_features.py ==========
    echo "Step 4: Merging surface features to PKL..."
    TMPDIR="$PROCESS_TMP_DIR" MASIF_TMP_DIR="$PROCESS_TMP_DIR" \
    python3 /es01/paratera/sce0413/czn/preprocess_pkl/preprocess/merge_surface_features.py \
        --pkl_file "$pkl_output" \
        --ligand_npz "$ligand_npz" \
        --pocket_npz "$pocket_npz" \
        --output_file "$output_dir/${name}_features_with_masif.pkl" < /dev/null
    if [ $? -ne 0 ]; then
        echo "Error in merge_surface_features.py for $name"
        return 1
    fi

    merged_pkl="$output_dir/${name}_features_with_masif.pkl"
    if [ ! -f "$merged_pkl" ]; then
        echo "Error: Merged PKL file not generated: $merged_pkl"
        return 1
    fi

    # ========== 清理中间文件 ==========
    echo "Cleaning up intermediate files..."
    rm -rf "$output_dir/tmp" 2>/dev/null || true
    rm -rf "$output_dir/surfaces" 2>/dev/null || true

    echo "Successfully processed: $name"
    return 0
}

export -f process_one_complex

<<<<<<< HEAD
LOG_DIR="/es01/paratera/sce0413/czn/log_update/$(basename "$CSV_FILE" .csv)"
mkdir -p "$LOG_DIR"

=======
>>>>>>> 4f9f0f7 (ppi_to_be_done)
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

success_count=$(grep -c "Successfully processed:" "$LOG_DIR/output_success.log" 2>/dev/null || echo "0")
fail_count=$(grep -c "Error" "$LOG_DIR/output_fail.log" 2>/dev/null || echo "0")

echo "=========================================="
echo "Processing finished!"
echo "Success: $success_count"
echo "Failed: $fail_count"
echo "Logs saved in: $LOG_DIR"
echo "=========================================="