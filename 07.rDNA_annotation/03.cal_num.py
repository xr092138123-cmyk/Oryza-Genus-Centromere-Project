import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="统计各种rDNA的数量")
parser.add_argument('--bed', type=str, required=True, help='Path to the BED file')
args = parser.parse_args()
bed = args.bed

def get_count(bed):
    """统计rDNA数量"""
    prefix = bed.rsplit('.',1)[0]
    species = bed.split('.')[0]
    temp = pd.read_csv(bed,header=None, sep = '\t')
    temp = temp.replace(r'Name=', '',regex=True)
    count = temp.groupby(0)[3].value_counts().unstack(fill_value=0)
    count['species'] = species
    output_name = f"{prefix}.num"
    output = count.loc[:,['18S_rRNA', '5_8S_rRNA', '28S_rRNA', '5S_rRNA', 'species']]
    output.to_csv(output_name, sep = '\t')
    print(f'*****已统计{bed}中各种rDNA的数量*****')

get_count(bed)