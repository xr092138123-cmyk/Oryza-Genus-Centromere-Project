#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import math
import sys

# 定义您指定的中性替代率 (substitutions/site/year)
NEUTRAL_RATE = 1.3e-8

def parse_fasta(filepath):
    """
    一个简单的FASTA文件解析器。
    
    参数:
        filepath (str): FASTA文件的路径。
    
    返回:
        dict: 一个字典，键是序列头部，值是序列字符串。
    """
    sequences = {}
    with open(filepath, 'r') as f:
        header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(current_seq)
                header = line[1:].strip() # strip() added to remove potential whitespace
                current_seq = []
            else:
                current_seq.append(line.upper())
        if header:
            sequences[header] = ''.join(current_seq)
    return sequences

def calculate_jukes_cantor_distance(seq1, seq2):
    """
    使用Jukes-Cantor (JC69)模型计算两条序列间的遗传距离 (K)。
    
    该模型假设所有碱基替换以相同速率发生。
    它会忽略包含'-'或'N'的位点。

    参数:
        seq1 (str): 第一条序列。
        seq2 (str): 第二条序列。

    返回:
        tuple: (遗传距离K, 原始差异p) 或 (None, None) 如果无法计算。
    """
    if len(seq1) != len(seq2):
        raise ValueError("错误：两条序列必须是对齐的，且长度相同。")

    diffs = 0
    compared_sites = 0
    valid_bases = 'ATGC'

    for base1, base2 in zip(seq1, seq2):
        # 只比较两个位点都是明确碱基的位置
        if base1 in valid_bases and base2 in valid_bases:
            compared_sites += 1
            if base1 != base2:
                diffs += 1
    
    if compared_sites == 0:
        return 0.0, 0.0 # 没有可比较的位点

    # p-distance: 观察到的差异比例
    p_distance = diffs / compared_sites

    # JC69模型的限制：如果p_distance >= 3/4，距离会趋于无穷大
    # 这表示序列差异太大（饱和），无法可靠估算
    if p_distance >= 0.75:
        print(
            f"警告：序列差异过大 (p-distance = {p_distance:.4f} >= 0.75)，"
            "无法使用Jukes-Cantor模型进行可靠计算。",
            file=sys.stderr
        )
        return None, p_distance
    
    # Jukes-Cantor (JC69) 公式: K = -3/4 * ln(1 - 4/3 * p)
    # 处理 p_distance = 0 的情况，避免 log(1) 计算浮点误差
    if p_distance == 0.0:
        return 0.0, 0.0
        
    k_distance = -0.75 * math.log(1.0 - (4.0 / 3.0) * p_distance)
    
    return k_distance, p_distance

def main():
    """主执行函数"""
    parser = argparse.ArgumentParser(
        description="计算NUMT的插入时间，并输出为TSV格式。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help="输入的FASTA文件路径。\n该文件必须包含两条对齐的序列：NUMT及其对应的线粒体序列。"
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help="输出结果的TSV文件路径。"
    )
    args = parser.parse_args()

    # 1. 解析FASTA文件
    try:
        sequences = parse_fasta(args.input)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{args.input}'", file=sys.stderr)
        sys.exit(1)

    if len(sequences) != 2:
        print(
            f"错误：输入文件 '{args.input}' 必须恰好包含两条序列，但检测到 {len(sequences)} 条。",
            file=sys.stderr
        )
        sys.exit(1)
        
    headers = list(sequences.keys())
    seq1, seq2 = list(sequences.values())

    # 2. 计算遗传距离 (K)
    try:
        k, _ = calculate_jukes_cantor_distance(seq1, seq2)
        if k is None:
            sys.exit(1) # 如果距离无法计算，则退出
    except ValueError as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    # 3. 计算插入时间 (T = K / r)
    insertion_time_years = k / NEUTRAL_RATE

    # 4. 准备TSV格式的输出行
    # 列: 序列1名称; 序列2名称; 中性替代率; 遗传距离模型; 校正后遗传距离; 插入时间（年）
    model_name = "Jukes-Cantor(JC69)"
    output_line = (
        f"{headers[0]}\t"
        f"{headers[1]}\t"
        f"{NEUTRAL_RATE:.2e}\t"
        f"{model_name}\t"
        f"{k:.6f}\t"
        f"{insertion_time_years:.2f}\n"
    )

    # 5. 将TSV行写入输出文件
    try:
        with open(args.output, 'w') as f:
            f.write(output_line)
        print(f"计算完成！TSV格式的结果已成功写入到 '{args.output}'。")
    except IOError as e:
        print(f"错误：无法写入输出文件 '{args.output}'。原因: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
