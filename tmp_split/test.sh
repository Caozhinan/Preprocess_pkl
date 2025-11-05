#!/bin/bash
# set -euo pipefail

OUT="merged_jobs.csv"

# 先写入 job_1.csv 的表头
head -n1 job_1.csv > "$OUT"

# 循环所有 job_*.csv 把数据部分拼接
for f in job_*.csv; do
    echo "Merging $f"
    tail -n +2 "$f" >> "$OUT"
done

echo "Merged file: $OUT"