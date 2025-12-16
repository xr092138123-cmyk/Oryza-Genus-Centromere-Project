import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import argparse
# 读取文件参数
parser = argparse.ArgumentParser(description="对rDNA的统计结果绘制热图")
parser.add_argument('-i', '--input', type=str, required=True, help='Path to the BED file')
args = parser.parse_args()
inputfile = args.input

def add_missing_chr(heat):
    """传入一个数据框（可能有缺失的染色体），返回一个有所有染色体的数据框"""
    Chr = ['Chr01', 'Chr02', 'Chr03', 'Chr04', 'Chr05', 'Chr06', 'Chr07', 'Chr08', 'Chr09', 'Chr10', 'Chr11', 'Chr12']
    # 检查 Chr 列表中的元素是否在 heat 的索引中
    missing_elements = [elem for elem in Chr if elem not in heat.index]
    
    # 如果有缺失的元素，则添加到 DataFrame 中
    if missing_elements:
        # 创建一个新 DataFrame，其中缺失的元素作为索引，所有值为 0
        new_rows = pd.DataFrame(0, index=missing_elements, columns=heat.columns)
        # 将新行添加到原始 DataFrame 中
        heat = pd.concat([heat, new_rows])
    heat = heat.sort_index()
    return heat

def make_heatmap(inputfile):
    prefix = inputfile.rsplit('.',1)[0]
    outputfile = f"{prefix}.heatmap.png"
    # 读取文件并处理
    temp = pd.read_csv(inputfile, sep = '\t')
    temp = temp.set_index('0')
    temp = temp.drop(columns = 'species')
    newtemp = add_missing_chr(temp)
    temp5s = newtemp.loc[:,['5S_rRNA']]
    temp48s = newtemp.drop(columns = '5S_rRNA')
    
    # 绘图
    # 自定义颜色映射表
    colors1 = ['#f4dbd9', '#cb364a']
    mycmap1 = LinearSegmentedColormap.from_list('custom_cmap', colors1)
    colors2 = ['#cfe3ee', '#5292c0']
    mycmap2 = LinearSegmentedColormap.from_list('custom_cmap', colors2)
    # 绘制热图
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(1,4)
    ax1 = fig.add_subplot(gs[0,:3])
    sns.heatmap(temp48s, ax=ax1, cmap=mycmap1, annot=True, fmt='.0f')
    ax2 = fig.add_subplot(gs[0, 3])
    sns.heatmap(temp5s, ax=ax2, cmap=mycmap2, annot=True, fmt='.0f')
    plt.savefig(outputfile, dpi = 300, bbox_inches='tight')

make_heatmap(inputfile)
print(f'*****已绘制{inputfile}的热图*****')
