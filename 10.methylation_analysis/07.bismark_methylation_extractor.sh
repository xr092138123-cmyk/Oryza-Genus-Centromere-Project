#!/bin/bash

# ===== 全局参数配置 =====
dedup_dir="/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/02_bismark.2/03_deduplicate_hap1"
output_base="/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/02_bismark.2/04_bismark_methylation_hap1"
index_base="/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/02_bismark.2/01_bismark_genome_preparation"
nip_index="/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/02_bismark.1/01_bismark_genome_preparation/NIP"
#QUEUE="c01"  # 队列名称
# 排除样本列表（使用前缀匹配）
exclude_samples=("R119" "R21CX-24" "R21CX-32" "R21CX-8" "NIP")

# ===== 自动检测BAM文件并提交作业 =====
for bam_file in ${dedup_dir}/*.deduplicated.bam; do
    # 提取基础文件名（不含路径和扩展名）
    base_name=$(basename "$bam_file" .deduplicated.bam)
    
    # 跳过排除样本
    skip=false
    for exclude in "${exclude_samples[@]}"; do
        if [[ "$base_name" == *"$exclude"* ]]; then
            echo "跳过排除样本: $base_name"
            skip=true
            break
        fi
    done
    if $skip; then continue; fi

    # 创建样本专属输出目录
    sample_outdir="${output_base}/${base_name}"
    mkdir -p "$sample_outdir"
    
    # 确定索引路径（特殊处理NIP样本）
    if [[ "$base_name" == Nip_* ]]; then
        genome_index="$nip_index"
    else
        # 提取样本前缀（如W3012_1 -> W3012）
        sample_prefix=$(echo "$base_name" | cut -d'_' -f1)
        genome_index="${index_base}/${sample_prefix}_hap1"
        
        # 验证索引目录存在
        if [ ! -d "$genome_index" ]; then
            echo "错误：索引目录不存在 $genome_index"
            exit 1
        fi
    fi

    # 生成作业提交脚本
    cat > "${sample_outdir}/submit_job.sh" << EOF
#!/bin/bash
#CSUB -J ${base_name}_methylation
#CSUB -q c01
#CSUB -o ${base_name}.o
#CSUB -e ${base_name}.e
#CSUB -n 8
#CSUB -R span[hosts=1]
#CSUB -cwd $sample_outdir

# 执行甲基化提取
bismark_methylation_extractor \
  -p \
  --bedGraph \
  --CX \
  --buffer_size 30G \
  --parallel 4 \
  --cytosine_report \
  --report \
  --counts \
  --comprehensive \
  --no_overlap \
  --genome_folder "$genome_index" \
  "$bam_file" \
  -o "$sample_outdir"

echo "完成样本: $base_name | 输出目录: $sample_outdir"
EOF

    # 提交作业并设置权限
    chmod u+x "${sample_outdir}/submit_job.sh"
    csub <"${sample_outdir}/submit_job.sh"
    echo "已提交作业: ${base_name}"
done

echo "所有作业提交完成! 结果将统一保存在: $output_base"
