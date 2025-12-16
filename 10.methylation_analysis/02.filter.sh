# ---------- 关键路径配置 ----------
INPUT_DIR="rice_methydata/WGBS.data/data"
OUTPUT_DIR="rice_methy_analyze/01_fastp"

# 提取所有唯一前缀（如HZ_1, R21CX-24_1等）
mapfile -t sample_ids < <(find "$INPUT_DIR" -name "*_clean_R1.fastq" -exec basename {} \; | sed 's/_clean_R1.fastq//')

# ---------- 并行处理循环 ----------
for id in "${sample_ids[@]}"; do
    # 构建文件路径（双引号包裹防特殊字符）
    r1_in="${INPUT_DIR}/${id}_clean_R1.fastq"
    r2_in="${INPUT_DIR}/${id}_clean_R2.fastq"
    r1_out="${OUTPUT_DIR}/${id}_clean_R1.fastq"
    r2_out="${OUTPUT_DIR}/${id}_clean_R2.fastq"
    report_html="${OUTPUT_DIR}/${id}_report.html"
    report_json="${OUTPUT_DIR}/${id}_report.json"

    # 检查文件是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "[ERROR] 文件缺失: $id (跳过处理)" >&2
        continue
    fi

    # 执行Fastp（后台并行）
    fastp -i "$r1_in" -I "$r2_in" \
          -o "$r1_out" -O "$r2_out" \
          -q 20 -l 50 \
          --detect_adapter_for_pe \
          -h "$report_html" -j "$report_json" \
          --thread 16 &  # 后台运行
done

# ---------- 等待所有任务完成 ----------
wait
echo "所有Fastp质控任务已完成！报告保存在: $OUTPUT_DIR"
...skipping...
    # 构建文件路径（双引号包裹防特殊字符）
    r1_in="${INPUT_DIR}/${id}_clean_R1.fastq"
    r2_in="${INPUT_DIR}/${id}_clean_R2.fastq"
    r1_out="${OUTPUT_DIR}/${id}_clean_R1.fastq"
    r2_out="${OUTPUT_DIR}/${id}_clean_R2.fastq"
    report_html="${OUTPUT_DIR}/${id}_report.html"
    report_json="${OUTPUT_DIR}/${id}_report.json"

    # 检查文件是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "[ERROR] 文件缺失: $id (跳过处理)" >&2
        continue
    fi

    # 执行Fastp（后台并行）
    fastp -i "$r1_in" -I "$r2_in" \
          -o "$r1_out" -O "$r2_out" \
          -q 20 -l 50 \
          --detect_adapter_for_pe \
          -h "$report_html" -j "$report_json" \
          --thread 16 &  # 后台运行
done

# ---------- 等待所有任务完成 ----------
wait
echo "所有Fastp质控任务已完成！报告保存在: $OUTPUT_DIR"
...skipping...
    # 构建文件路径（双引号包裹防特殊字符）
    r1_in="${INPUT_DIR}/${id}_clean_R1.fastq"
    r2_in="${INPUT_DIR}/${id}_clean_R2.fastq"
    r1_out="${OUTPUT_DIR}/${id}_clean_R1.fastq"
    r2_out="${OUTPUT_DIR}/${id}_clean_R2.fastq"
    report_html="${OUTPUT_DIR}/${id}_report.html"
    report_json="${OUTPUT_DIR}/${id}_report.json"

    # 检查文件是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "[ERROR] 文件缺失: $id (跳过处理)" >&2
        continue
    fi

    # 执行Fastp（后台并行）
    fastp -i "$r1_in" -I "$r2_in" \
          -o "$r1_out" -O "$r2_out" \
          -q 20 -l 50 \
          --detect_adapter_for_pe \
          -h "$report_html" -j "$report_json" \
          --thread 16 &  # 后台运行
done

# ---------- 等待所有任务完成 ----------
wait
echo "所有Fastp质控任务已完成！报告保存在: $OUTPUT_DIR"
...skipping...
    # 构建文件路径（双引号包裹防特殊字符）
    r1_in="${INPUT_DIR}/${id}_clean_R1.fastq"
    r2_in="${INPUT_DIR}/${id}_clean_R2.fastq"
    r1_out="${OUTPUT_DIR}/${id}_clean_R1.fastq"
    r2_out="${OUTPUT_DIR}/${id}_clean_R2.fastq"
    report_html="${OUTPUT_DIR}/${id}_report.html"
    report_json="${OUTPUT_DIR}/${id}_report.json"

    # 检查文件是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "[ERROR] 文件缺失: $id (跳过处理)" >&2
        continue
    fi

    # 执行Fastp（后台并行）
    fastp -i "$r1_in" -I "$r2_in" \
          -o "$r1_out" -O "$r2_out" \
          -q 20 -l 50 \
          --detect_adapter_for_pe \
          -h "$report_html" -j "$report_json" \
          --thread 16 &  # 后台运行
done

# ---------- 等待所有任务完成 ----------
wait
echo "所有Fastp质控任务已完成！报告保存在: $OUTPUT_DIR"
...skipping...
    # 构建文件路径（双引号包裹防特殊字符）
    r1_in="${INPUT_DIR}/${id}_clean_R1.fastq"
    r2_in="${INPUT_DIR}/${id}_clean_R2.fastq"
    r1_out="${OUTPUT_DIR}/${id}_clean_R1.fastq"
    r2_out="${OUTPUT_DIR}/${id}_clean_R2.fastq"
    report_html="${OUTPUT_DIR}/${id}_report.html"
    report_json="${OUTPUT_DIR}/${id}_report.json"

    # 检查文件是否存在
    if [[ ! -f "$r1_in" || ! -f "$r2_in" ]]; then
        echo "[ERROR] 文件缺失: $id (跳过处理)" >&2
        continue
    fi

    # 执行Fastp（后台并行）
    fastp -i "$r1_in" -I "$r2_in" \
          -o "$r1_out" -O "$r2_out" \
          -q 20 -l 50 \
          --detect_adapter_for_pe \
          -h "$report_html" -j "$report_json" \
          --thread 16 &  # 后台运行
done

# ---------- 等待所有任务完成 ----------
wait
echo "所有Fastp质控任务已完成！报告保存在: $OUTPUT_DIR"