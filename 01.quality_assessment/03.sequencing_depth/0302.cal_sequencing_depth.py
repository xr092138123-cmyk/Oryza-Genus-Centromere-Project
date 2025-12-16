import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="对文件进行排序")
parser.add_argument('--ont', type=str, required=True, help='测序数据文件')
args = parser.parse_args()

ont = args.ont

# 输入文件类似于AA_Oruf.ont.sorted.bam


prefix = ont.split('.')[0]
# 建立索引
command5 = f"samtools index {prefix}.ont.sorted.bam"
os.system(command5)
# 计算全基因组范围内的覆盖度
command6 = f"mosdepth -t 64 {prefix}.hifi {prefix}.ont.sorted.bam"
os.system(command6)