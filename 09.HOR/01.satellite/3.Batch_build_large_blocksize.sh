#!/bin/bash
#CSUB -J HOR_blocks
#CSUB -q c02
#CSUB -o /share/org/YZWL/yzwl_shahd/WangJR/logs/%J.out
#CSUB -e /share/org/YZWL/yzwl_shahd/WangJR/logs/%J.err
#CSUB -n 1
#CSUB -cwd /share/org/YZWL/yzwl_shahd/WangJR/HORpairscore_from_Pairswise_Identity/
#CSUB -R span[hosts=1]

# 设置环境
module purge
module load python/3.9  # 根据实际调整版本

# 执行Python脚本（使用8个核心进行并行处理）
python /share/org/YZWL/yzwl_shahd/WangJR/HORpairscore_from_Pairswise_Identity/Batch_build_large_blocksize.py \
  --indir /share/org/YZWL/yzwl_shahd/WangJR/HORpairscore_from_Pairswise_Identity \
  --outdir /share/org/YZWL/yzwl_shahd/WangJR/large3 \
  --glob   "*.hor_pairs.tsv" \