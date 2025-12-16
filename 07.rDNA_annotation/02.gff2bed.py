import os
import pandas as pd
import argparse

# 再使用之前，需要在命令行使用bedops软件的gff2bed命令，将原始输出的gff文件转换为bed文件

parser = argparse.ArgumentParser(description="将bed文件转换为标准的bed文件")
parser.add_argument('--bed', type=str, required=True, help='Path to the BED file')
args = parser.parse_args()
bed = args.bed

def get_standard_bed(bed):
    """输入一个gff文件转换后的bed文件，得到一个标准的6列的bed文件"""
    material = bed.split('.')[0]
    temp = pd.read_csv(bed, header=None, sep='\t|;', engine='python')
    output = temp.loc[:, [0, 1, 2, 9, 4, 5]]
    output[9] = output[9].replace(r'Name=', '', regex=True)
    #output[9] = output[9].replace(r'rRNA', 'rDNA', regex=True)
    #output[9] = output[9] + '_' + (output[9].index + 1).astype(str)
    #output[9] = material + '-' + output[9]
    output_name = f"{bed}.txt"
    output.to_csv(output_name, header=False, sep = '\t', index=False)
    print(f"*****已将{bed}文件转换为标准的bed文件*****")

get_standard_bed(bed)
