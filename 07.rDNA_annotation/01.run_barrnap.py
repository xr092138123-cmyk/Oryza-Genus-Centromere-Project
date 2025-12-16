import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="对基因组进行rDNA注释")
parser.add_argument('-fa', '--fasta', type=str, required=True, help='基因组文件')
args = parser.parse_args()

def barrnap(fa):
    """输入fasta序列文件，输出rDNA的序列以及注释文件"""
    prefix = fa.rsplit('.',1)[0]
    command1 = f"pybarrnap -k euk {fa} --accurate > {prefix}.rDNA.gff"
    os.system(command1)
    print(f"*****已对{fa}的rDNA进行注释*****")

barrnap(args.fasta)
