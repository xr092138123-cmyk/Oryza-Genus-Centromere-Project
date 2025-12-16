import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="将测序数据比对到基因组上")
parser.add_argument('-d', '--data', type=str, required=True, help='测序数据类型ont或hifi')
parser.add_argument('-s', '--seq', type=str, required=True, help='测序数据类型ont或hifi')
parser.add_argument('-g', '--genome', type=str, required=True, help='基因组数据')
args = parser.parse_args()

def hifi_mapping_coverage(hifi, genome):
    """hifi数据比对"""
    prefix = hifi.split('.')[0]
    # 将原始的bam格式转换为fastq格式
    command1 = f"bedtools bamtofastq -i {hifi} -fq {prefix}.hifi.fastq"
    os.system(command1)
    # 比对
    command2 = f"minimap2 -ax map-pb {genome} {prefix}.hifi.fastq -t 64 > {prefix}.hifi.sam"
    os.system(command2)
    # 将SAM文件转换为BAM文件
    command3 = f"samtools view -@ 64 -Sb {prefix}.hifi.sam > {prefix}.hifi.bam"
    os.system(command3)
    # 排序
    command4 = f"samtools sort -@ 64 -o {prefix}.hifi.sorted.bam {prefix}.hifi.bam"
    os.system(command4)

def ont_mapping_coverage(ont, genome):
    """ont数据比对"""
    prefix = ont.split('.')[0]
    # 比对
    command1 = f"minimap2 --secondary=no -ax map-ont {genome} {ont} -t 64 > {prefix}.ont.sam"
    #os.system(command1)
    # 将sam文件转换为bam文件
    command2 = f"samtools view -@ 64 -Sb {prefix}.ont.sam > {prefix}.ont.bam"
    os.system(command2)
    # 排序
    command3 = f"samtools sort -@ 64 -o {prefix}.ont.sorted.bam {prefix}.ont.bam"
    os.system(command3)


# 执行
if args.data == 'hifi':
    hifi_mapping_coverage(args.seq, args.genome)
elif args.data == 'ont':
    ont_mapping_coverage(args.seq, args.genome)