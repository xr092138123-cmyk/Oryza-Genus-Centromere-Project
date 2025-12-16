import pandas as pd
import argparse
import sys

def filter_blast_results(blast_file, output_file, identity_threshold=85.0, evalue_threshold=1e-5):
    """
    使用pandas筛选BLAST outfmt 6的结果。

    Args:
        blast_file (str): 输入的BLAST文件名 (outfmt 6)。
        output_file (str): 输出的筛选后文件名。
        identity_threshold (float): 百分比同一性的最小阈值。
        evalue_threshold (float): E-value的最大阈值。
    """
    # 定义BLAST outfmt 6的标准列名，方便后续引用
    column_names = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]

    try:
        print(f"正在读取 BLAST 文件: {blast_file} ...")
        # 使用pandas读取制表符分隔的文件，并指定列名
        df = pd.read_csv(
            blast_file, 
            sep='\t',         # 分隔符是制表符
            header=None,      # 文件没有标题行
            names=column_names # 指定我们定义的列名
        )
        
        num_original_records = len(df)
        print(f"-> 读取完成，共找到 {num_original_records} 条比对记录。")

    except FileNotFoundError:
        print(f"错误: 输入文件 '{blast_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"错误: 输入文件 '{blast_file}' 为空，无法处理。", file=sys.stderr)
        sys.exit(1)

    # --- 核心筛选逻辑 ---
    print(f"\n开始筛选...")
    print(f"筛选条件1: identity (pident) > {identity_threshold}")
    print(f"筛选条件2: E-value (evalue) < {evalue_threshold}")

    # 使用布尔索引进行筛选，'&' 表示 '与' (AND)
    filtered_df = df[
        (df['pident'] > identity_threshold) & 
        (df['evalue'] < evalue_threshold)
    ]
    
    num_filtered_records = len(filtered_df)
    print(f"-> 筛选完成，剩下 {num_filtered_records} 条记录。")

    # --- 保存结果 ---
    try:
        print(f"\n正在将筛选结果保存到: {output_file} ...")
        # 将筛选后的DataFrame保存为新的制表符分隔文件
        # index=False: 不保存pandas的行索引
        # header=False: 不写入列名，保持与原始outfmt 6格式一致
        filtered_df.to_csv(output_file, sep='\t', index=False, header=False)
        
        print("-" * 30)
        print("处理成功！")
        print(f"原始记录数: {num_original_records}")
        print(f"筛选后记录数: {num_filtered_records}")
        print(f"结果已保存至: {output_file}")
        print("-" * 30)
        
        # 打印前5行结果作为预览
        print("筛选结果预览 (前5行):")
        # to_string() 可以更好地格式化输出
        print(filtered_df.head().to_string(index=False, header=False))

    except Exception as e:
        print(f"错误: 保存文件时发生错误: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    主函数，用于解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="使用pandas筛选BLAST outfmt 6的结果文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--blast_file", 
        help="输入的BLAST结果文件 (outfmt 6)。\n例如: blast.txt"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="输出的筛选后文件名。\n例如: filtered_blast.txt"
    )
    parser.add_argument(
        "--identity", 
        type=float, 
        default=85.0, 
        help="Identity的最小阈值 (默认: 85.0)"
    )
    parser.add_argument(
        "-e", "--evalue", 
        type=float, 
        default=1e-5, 
        help="E-value的最大阈值 (默认: 1e-5)"
    )
    
    args = parser.parse_args()
    
    filter_blast_results(args.blast_file, args.output_file, args.identity, args.evalue)

if __name__ == "__main__":
    main()
