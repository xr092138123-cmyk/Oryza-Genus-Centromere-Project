#!/usr/bin/env python3

import argparse
import sys
from itertools import combinations

# 尝试导入Levenshtein库，如果失败则给出提示
try:
    import Levenshtein
except ImportError:
    print("错误：缺少 'Levenshtein' 库。", file=sys.stderr)
    print("请使用 'pip install python-Levenshtein' 命令进行安装。", file=sys.stderr)
    sys.exit(1)

def parse_fasta(file_path):
    """
    解析FASTA文件，返回一个序列名称到序列的字典。
    
    Args:
        file_path (str): FASTA文件的路径。
        
    Returns:
        dict: 一个字典，键是序列名，值是DNA序列。
    """
    sequences = {}
    current_header = None
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    current_header = line[1:]
                    sequences[current_header] = []
                elif current_header:
                    sequences[current_header].append(line)

        for header, seq_list in sequences.items():
            sequences[header] = "".join(seq_list).upper()
            
        return sequences
    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{file_path}'", file=sys.stderr)
        sys.exit(1)

def calculate_circular_min_distance(seq1, seq2):
    """
    计算两条序列间的最小循环编辑距离。
    
    将seq2进行所有可能的循环位移，分别与固定的seq1计算编辑距离，
    并返回其中的最小值。即使序列长度不同，此操作也会执行。
    
    Args:
        seq1 (str): 第一条序列（固定的参考序列）。
        seq2 (str): 第二条序列（将被循环位移的序列）。
        
    Returns:
        int: 最小的编辑距离。
    """
    min_distance = float('inf')
    temp_seq2 = seq2
    
    if not temp_seq2:
        return Levenshtein.distance(seq1, seq2)

    for _ in range(len(seq2)):
        dist = Levenshtein.distance(seq1, temp_seq2)
        if dist < min_distance:
            min_distance = dist
        if min_distance == 0:
            break
        temp_seq2 = temp_seq2[1:] + temp_seq2[0]
        
    return min_distance

def main():
    """主执行函数"""
    parser = argparse.ArgumentParser(
        description="计算FASTA文件中序列对之间的最小循环编辑距离（支持不同长度序列，流式写入以节省内存）。"
    )
    parser.add_argument(
        "--input", 
        required=True, 
        help="输入的FASTA文件路径。"
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="输出结果的TSV文件路径。"
    )
    
    args = parser.parse_args()
    
    # 1. 解析FASTA文件
    print(f" 正在读取并解析文件: {args.input}")
    sequences = parse_fasta(args.input)
    if not sequences or len(sequences) < 2:
        print("错误：输入文件必须包含至少两条有效的序列才能进行比较。", file=sys.stderr)
        sys.exit(1)
    print(f" 成功读取 {len(sequences)} 条序列。")
    
    print(" 正在计算并流式写入结果到: {}...".format(args.output))
    
    try:
        # 2. 在循环开始前打开输出文件
        with open(args.output, 'w') as f_out:
            # 3. 首先写入表头
            f_out.write("Sequence_A\tSequence_B\tMin_Edit_Distance\n")
            
            sequence_names = list(sequences.keys())
            pair_iterator = combinations(sequence_names, 2)
            
            # 使用一个计数器来估算总数并显示进度
            # 对于非常大的N，直接计算combinations的长度可能会很慢或消耗内存
            # 所以我们用一个简单的方法来估算
            num_seqs = len(sequence_names)
            total_pairs = num_seqs * (num_seqs - 1) // 2
            
            # 4. 进入计算循环
            for i, (name1, name2) in enumerate(pair_iterator):
                seq1 = sequences[name1]
                seq2 = sequences[name2]
                
                # 计算核心逻辑
                min_dist = calculate_circular_min_distance(seq1, seq2)
                
                # 5. 每计算完一对，就立即写入文件
                f_out.write(f"{name1}\t{name2}\t{min_dist}\n")
                
                # 打印进度
                if (i + 1) % 100 == 0 or (i + 1) == total_pairs:
                    progress = (i + 1) / total_pairs * 100
                    print(f"  ...已处理 {i+1}/{total_pairs} 对 ({progress:.1f}%)", end='\r')

    except IOError as e:
        print(f"\n错误：无法写入输出文件 '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)
        
    print("\n\n 所有任务完成！")

if __name__ == "__main__":
    main()