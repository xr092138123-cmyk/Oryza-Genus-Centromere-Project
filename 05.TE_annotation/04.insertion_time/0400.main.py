import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="计算TE的插入时间")
parser.add_argument('--bed', type=str, required=True, help='HiTE记录所有TE的bed文件')
parser.add_argument('--genome', type=str, required=True, help='基因组文件')
parser.add_argument('--te', type=str, required=True, help='转座子的文件')
args = parser.parse_args()

def main(bed, genome, te):
    prefix = bed.split('.')[0]
    # 将每对LTR进行拆分，每队LTR保存为一个bed文件，该bed文件中有lLTR和rLTR
    command_paired_LTR = f"python 01.split_paired_LTR.py --input {bed} --output_dir {prefix}.paired_LTR"
    os.system(command_paired_LTR)
    command_getfasta = f"for file in {prefix}.paired_LTR/*.bed; do bedtools getfasta -fi {genome} -fo ${{file/.bed/.fa}} -bed ${{file}} -s -nameOnly; done"
    os.system(command_getfasta)
    command_mafft = f"for file in {prefix}.paired_LTR/*.fa; do mafft --auto ${{file}} > ${{file/.fa/.mafft}}; done"
    os.system(command_mafft)

    # 计算插入时间，输出文件的每一列分别为：序列1名称；序列2名称；中性替代率；遗传距离模型；校正后遗传距离 (K)；插入时间（年）
    command_time = f"for file in {prefix}.paired_LTR/*.mafft; do python 02.insertion_time.py --input ${{file}} --output ${{file/.mafft/.time}}; done"
    os.system(command_time)

    # 给插入时间添加一列，该列为repeat region的名称
    command_time_addinfo = f"for file in {prefix}.paired_LTR/*.time; do python 03.add_info.py --input ${{file}}; done"
    os.system(command_time_addinfo)
    # 合并插入时间
    command_cat = f"cat {prefix}.paired_LTR/*.time > all.{prefix}.time"
    os.system(command_cat)

    # 将插入时间比对回intactLTR的bed文件中去
    command_map_time2bed = f"python 04.map_time2bed.py --bed {te} --time all.{prefix}.time --output {te}.withtime"
    os.system(command_map_time2bed)


main(args.bed, args.genome, args.te)

