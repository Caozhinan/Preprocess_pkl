#!/bin/bash
#SBATCH -J merge_pkl_parallel
#SBATCH -p G1Part_sce
#SBATCH -N 1
#SBATCH -n 1              # 使用 12 核
#SBATCH -t 48:00:00
#SBATCH -o del_dirs_all.out
#SBATCH -e del_dirs_all.err
# set -euo pipefail

LIST_FILE="/es01/paratera/sce0413/czn/preprocess_pkl/tmp_split/abs_dirs.txt"
WORK_DIR="/es01/paratera/sce0413/czn/preprocess_pkl/tmp_split"
SPLIT_PREFIX="$WORK_DIR/part_"
NUM_PARTS=18

if [ ! -f "$LIST_FILE" ]; then
    echo "Error: $LIST_FILE not found"
    exit 1
fi

echo "Splitting $LIST_FILE into $NUM_PARTS parts..."

# 获取总行数
TOTAL_LINES=$(wc -l < "$LIST_FILE")
LINES_PER_PART=$(( (TOTAL_LINES + NUM_PARTS - 1) / NUM_PARTS ))

# 按行数拆分
split -l $LINES_PER_PART "$LIST_FILE" "$SPLIT_PREFIX"

echo "Split done. Generated parts:"
ls ${SPLIT_PREFIX}*

# 删除脚本，写到 $WORK_DIR/delete_dirs.sh
DEL_SCRIPT="$WORK_DIR/delete_dirs.sh"
cat > "$DEL_SCRIPT" <<'EOS'
#!/bin/bash
set -euo pipefail
LIST_FILE="$1"
echo "Processing $LIST_FILE ..."
while IFS= read -r dir; do
    if [ -d "$dir" ]; then
        echo "Deleting $dir"
        rm -rf "$dir"
    else
        echo "Skip (not found): $dir"
    fi
done < "$LIST_FILE"
echo "Finished $LIST_FILE"
EOS
chmod +x "$DEL_SCRIPT"

# 遍历每个 part_* 提交 sbatch
for f in ${SPLIT_PREFIX}*; do
    sbatch <<EOF
#!/bin/bash
#SBATCH -J del_$(basename $f)
#SBATCH -o ${f}.out
#SBATCH -e ${f}.err
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0:30:00

bash "$DEL_SCRIPT" "$f"
EOF
done

echo "Submitted deletion jobs for all parts."