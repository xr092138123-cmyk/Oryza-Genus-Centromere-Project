#!/bin/bash
#CSUB -J iqtree
#CSUB -q c02
#CSUB -o /share/org/YZWL/yzwl_shahd/WangJR/logs/%J.out
#CSUB -e /share/org/YZWL/yzwl_shahd/WangJR/logs/%J.err
#CSUB -n 64
#CSUB -cwd /share/org/YZWL/yzwl_shahd/WangJR/HOR/Script
#CSUB -R span[hosts=1]

# 设置环境
module purge
module load python/3.9  # 根据实际调整版本

# 创建必要的目录
mkdir -p /share/org/YZWL/yzwl_shahd/WangJR/HOR/HORpairscore_from_Pairswise_Identity_maxblock50
mkdir -p /share/org/YZWL/yzwl_shahd/WangJR/logs

# 执行Python脚本
python Calculate_HORpairscore_from_Pairswise_Identity.py \
    --fa_glob "/share/org/YZWL/yzwl_shahd/wangdan/Oryza_Centromere/moddotplot/all_batch/06.repeat_top1_fa_new/*.fa" \
    --pairs_dir "/share/org/YZWL/yzwl_shahd/WangJR/HOR/Pairwise_Identity_from_moddot/" \
    --pairs_suffix ".pairwise.tsv" \
    --outdir /share/org/YZWL/yzwl_shahd/WangJR/HOR/HORpairscore_from_Pairswise_Identity_maxblock50 \
    --id_threshold 0.96 \
    --scale auto \
    --min_block 3 \
    --max_block 50 \
    --max_shift 0

# 检查执行状态
if [ $? -eq 0 ]; then
    echo "作业执行成功完成"
else
    echo "作业执行失败"
    exit 1
fi