import os
import argparse
import pandas as pd
import math
parser = argparse.ArgumentParser(description="QV值评估")
parser.add_argument('-g', '--genome', type=str, required=True, help='基因组数据')
args = parser.parse_args()
genome = args.genome

prefix = genome.split('.')[0]
## 求最佳k-mer数
command_seqkit = f"seqkit stats -a {genome} > temp.txt"
os.system(command_seqkit)
command_sed = fr"sed -i 's/ \+/\t/g' temp.txt"
os.system(command_sed)
temp = pd.read_csv('temp.txt', sep = '\t')
genome_len = temp.at[0, 'sum_len']
genome_len = genome_len.replace(",", "")
command_best_k = f"sh /share/org/YZWL/yzbsl_chengch/aaa_sand/merqury-master/best_k.sh {genome_len} > temp.txt"
os.system(command_best_k)
# 此时temp.txt中储存有最佳k-mer数，在第三行
# 打开文件并读取所有行
with open('temp.txt', 'r') as file:
    lines = file.readlines()
# 获取第三行（索引从 0 开始）
third_line = lines[2].strip()
third_line = float(third_line)
best_kmer = math.ceil(third_line)
print(f"*****该基因组{genome}的最佳k-mer为{best_kmer}")
os.remove('temp.txt')

# 创建k-mer数据集
command_meryl_R1 = f"meryl k={best_kmer} count output {prefix}.R1.meryl {prefix}.R1.fastq.gz"
command_meryl_R2 = f"meryl k={best_kmer} count output {prefix}.R2.meryl {prefix}.R2.fastq.gz"
os.system(command_meryl_R1)
os.system(command_meryl_R2)
command_meryl_merge = f"meryl union-sum output {prefix}.meryl {prefix}.R1.meryl {prefix}.R2.meryl"
os.system(command_meryl_merge)

# 运行merqury
command_merqury = f"sh /share/org/YZWL/yzbsl_chengch/aaa_sand/merqury-master/merqury.sh {prefix}.meryl {genome} {prefix}"
os.system(command_merqury)
