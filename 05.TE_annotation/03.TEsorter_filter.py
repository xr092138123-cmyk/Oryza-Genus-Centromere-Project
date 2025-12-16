import os
import time
import argparse
import pandas as pd
import re
import sys
parser = argparse.ArgumentParser(description="处理HiTE输出的LTR的内容，为intactLTR统一命名")
parser.add_argument('--fa', type=str, required=True, help='基因组序列文件')
parser.add_argument('--bed', type=str, required=True, help='intactLTR的注释文件，bed文件')
args = parser.parse_args()


def intactLTR_rename(bed, fa):
    """为intactLTR重新命名"""
    prefix = bed.split('.')[0]
    intact_fa = f"{prefix}.intactLTR.fa"
    command_getfasta = f"bedtools getfasta -fi {fa} -bed {bed} -fo temp.fa -s -nameOnly"
    os.system(command_getfasta)
    print(f"*****已提取出{intact_fa}的序列*****")
    clean_fasta_headers('temp.fa', intact_fa)
    os.remove('temp.fa')
    # 使用tesorter为intactLTR重命名
    tesorter_rename(intact_fa, bed)

def clean_fasta_headers(input_file, output_file):
    """
    读取一个FASTA文件，清理序列名称，然后写入一个新的FASTA文件。
    清理规则: 移除序列名称末尾的 `(?)`, `(+)`, 或 `(-)` 部分。

    Args:
        input_file (str): 输入的FASTA文件路径。
        output_file (str): 清理后输出的FASTA文件路径。
    """
    # 这个正则表达式匹配行尾的以下模式:
    # \s*     - 零个或多个空格
    # \(      - 一个字面的左括号
    # [?+\-]  - 括号内的一个字符，可以是 '?', '+', 或 '-'
    # \)      - 一个字面的右括号
    # $       - 表示匹配必须在行尾
    # 我们使用 re.compile() 来预编译正则表达式，这样在循环中使用时会更高效。
    pattern_to_remove = re.compile(r'\s*\([?+\-]\)$')

    print(f"开始处理文件: {input_file}")
    lines_processed = 0
    headers_cleaned = 0

    try:
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            for line in f_in:
                lines_processed += 1
                # 检查是否是序列名称行（以'>'开头）
                if line.startswith('>'):
                    # 移除行尾的换行符，以便进行处理
                    header = line.strip()
                    # 使用正则表达式的 sub() 方法替换匹配到的模式为空字符串
                    cleaned_header = pattern_to_remove.sub('', header)
                    
                    if header != cleaned_header:
                        headers_cleaned += 1
                    
                    # 将清理后的序列名称写回文件，并加上换行符
                    f_out.write(cleaned_header + '\n')
                else:
                    # 如果不是序列名称行（即序列本身），则直接写入
                    f_out.write(line)

    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理过程中发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n处理完成！")
    print(f"总共处理了 {lines_processed} 行。")
    print(f"清理了 {headers_cleaned} 个序列名称。")
    print(f"结果已保存到: {output_file}")

def tesorter_rename(input_fa, input_bed):
    """输入序列文件，及其对应的bed文件，将bed文件中的name替换为tesorter中注释到的name"""
    prefix = input_bed.rsplit('.',1)[0]
    # TEsorter注释, 获得*.cls.tsv注释文件
    command_tesorter = f"TEsorter {input_fa} -p 16 -pre {prefix} -db rexdb-plant"
    os.system(command_tesorter)
    # 将*.cls.tsv转换为ID~tsName的形式
    cls_tsv = f"{prefix}.cls.tsv"
    temp = pd.read_csv(cls_tsv, sep = '\t')
    temp['tsName'] = temp['Order'] + '_' + temp['Superfamily'] + '_' + temp['Clade']
    id_tsName = temp.loc[:, ['#TE','tsName']]
    id_tsName.columns = ['ID', 'tsName']
    id_tsName_map = f"{prefix}.tesorterLTR.txt"
    id_tsName.to_csv(id_tsName_map, index=False, sep = '\t')
    # 调用函数，进行转换
    id_tsName_bed = f"{prefix}.tesorterLTR.bed"
    ID_tsName_convert(input_bed, id_tsName_map, id_tsName_bed)
    print('-'*60)

def ID_tsName_convert(bed_path, map_path, output_path):
    """
    使用 pandas 处理 BED 文件，根据映射文件替换和重排 ID 列。

    Args:
        bed_path (str): 输入的 BED 文件路径。
        map_path (str): ID 到 tsName 的映射文件路径。
        output_path (str): 输出的修改后的文件路径。
    """
    try:
        # 1. 读取 BED 文件为 DataFrame，手动指定列名
        print(f"正在读取 BED 文件: {bed_path}")
        bed_columns = ['chrom', 'start', 'end', 'id', 'score', 'strand', 'anno']
        bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=bed_columns)
        
        # 2. 读取映射文件为 DataFrame，使用正则表达式匹配分隔符
        print(f"正在读取映射文件: {map_path}")
        map_df = pd.read_csv(map_path, sep='\t')

    except FileNotFoundError as e:
        print(f"错误: 文件未找到 - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

    print("文件读取成功，开始合并数据...")
    # 3. 合并两个 DataFrame
    # 使用左连接（left join）保留所有 bed_df 的行
    # left_on='id' 指的是 bed_df 中的 'id' 列
    # right_on='ID' 指的是 map_df 中的 'ID' 列 (根据文件中的表头)
    merged_df = pd.merge(bed_df, map_df, left_on='id', right_on='ID', how='left')

    # 4. 处理未匹配的 ID
    # 统计有多少 ID 没有在映射文件中找到
    unmapped_count = merged_df['tsName'].isnull().sum()
    if unmapped_count > 0:
        print(f"警告: 有 {unmapped_count} 个 ID 未能在映射文件中找到对应项。")
        # 如果 tsName 为空 (NaN)，则用原始 id 的值填充(改为用LTR_unknow_unknow来填充)
        merged_df['tsName'].fillna('none', inplace=True)

    # 5. 选择并重排最终的列
    # 新的第4列是 tsName，新的第5列是ID
    # 原始的 strand 和 map_df 的 ID 列被丢弃
    final_df = merged_df[['chrom', 'start', 'end', 'tsName', 'id', 'strand', 'anno']]

    # 6. 将结果写入文件
    print(f"正在将结果写入文件: {output_path}")
    try:
        final_df.to_csv(
            output_path, 
            sep='\t',        # 使用制表符作为分隔符
            header=False,    # 不写入列名（表头）
            index=False      # 不写入 DataFrame 的行索引
        )
    except Exception as e:
        print(f"写入文件时出错: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n处理完成！")
    print(f"成功处理 {len(bed_df)} 行 BED 数据。")
    print(f"结果已保存到: {output_path}")



intactLTR_rename(args.bed, args.fa)


