import pandas as pd
import argparse
import sys

def merge_bed_file(input_file, output_file):
    """
    读取、排序并合并BED文件中有重叠的区间。
    合并过程考虑链方向。

    Args:
        input_file (str): 输入的BED文件名。
        output_file (str): 输出的合并后的BED文件名。
    """
    # 定义BED文件的列名
    bed_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']

    try:
        print(f"正在读取 BED 文件: {input_file} ...")
        # 读取制表符分隔的BED文件
        df = pd.read_csv(
            input_file, 
            sep='\t', 
            header=None, 
            names=bed_columns,
            # 确保start和end被正确读取为整数
            dtype={'chrom': str, 'start': int, 'end': int, 'name': str, 'score': float, 'strand': str}
        )
        print(f"-> 读取完成，共找到 {len(df)} 条记录。")
    except FileNotFoundError:
        print(f"错误: 输入文件 '{input_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 关键步骤 1: 排序 ---
    # 必须先按 chrom, then strand, then start 排序
    print("正在按 'chrom', 'strand', 'start' 排序...")
    df_sorted = df.sort_values(by=['chrom', 'strand', 'start'])

    # --- 关键步骤 2: 合并重叠区间 ---
    print("正在合并重叠区间...")
    merged_records = []
    
    # 按 'chrom' 和 'strand' 分组处理
    for (chrom, strand), group in df_sorted.groupby(['chrom', 'strand']):
        
        # 获取当前组的第一个区间作为合并的起点
        if group.empty:
            continue
        
        # 使用 .iloc[0] 获取第一行的数据
        current_start = group.iloc[0]['start']
        current_end = group.iloc[0]['end']
        current_names = [group.iloc[0]['name']]
        current_scores = [group.iloc[0]['score']]

        # 从第二个区间开始遍历
        for i in range(1, len(group)):
            next_row = group.iloc[i]
            next_start = next_row['start']
            
            # 检查是否有重叠 (下一个区间的起始点在前一个区间的结束点之前)
            if next_start < current_end:
                # 有重叠，扩展当前合并区间
                current_end = max(current_end, next_row['end'])
                current_names.append(next_row['name'])
                current_scores.append(next_row['score'])
            else:
                # 没有重叠，完成上一个区间的合并，并将其写入结果列表
                merged_records.append({
                    'chrom': chrom,
                    'start': current_start,
                    'end': current_end,
                    'name': ','.join(map(str, sorted(list(set(current_names))))), # 合并name并去重
                    'score': max(current_scores), # 取最大score
                    'strand': strand
                })
                
                # 开始一个新的合并区间
                current_start = next_row['start']
                current_end = next_row['end']
                current_names = [next_row['name']]
                current_scores = [next_row['score']]

        # 不要忘记将最后一个合并的区间写入结果列表
        merged_records.append({
            'chrom': chrom,
            'start': current_start,
            'end': current_end,
            'name': ','.join(map(str, sorted(list(set(current_names))))),
            'score': max(current_scores),
            'strand': strand
        })

    # --- 步骤 3: 保存结果 ---
    if not merged_records:
        print("没有记录可以合并或输出。")
        # 创建一个空文件
        open(output_file, 'w').close()
        return
        
    print(f"-> 合并完成，生成 {len(merged_records)} 条新记录。")
    # 将结果列表转换为 DataFrame
    df_merged = pd.DataFrame(merged_records)
    
    # 按照最终结果再次排序（可选，但推荐）
    df_merged_sorted = df_merged.sort_values(by=['chrom', 'start'])
    
    print(f"正在将结果保存到: {output_file} ...")
    df_merged_sorted.to_csv(output_file, sep='\t', header=False, index=False)
    
    print("\n处理成功！")

def main():
    """
    主函数，用于解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="排序并合并BED文件中的重叠区间，同时考虑链方向。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--input_file", 
        help="输入的BED文件名。\n例如: input.bed"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="输出的合并后的BED文件名。\n例如: merged.bed"
    )
    
    args = parser.parse_args()
    merge_bed_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
