#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse

def parse_bismark_filename(filename):
    """
    根据用户定义的简单规则解析文件名。
    规则: 以下划线为分隔，第一个部分是品种名，第二个部分是重复数。

    Args:
        filename (str): 不含路径的文件名。

    Returns:
        tuple: (品种名, 重复数) 或 (None, None) 如果格式不匹配。
    """
    try:
        parts = filename.split('_')
        if len(parts) >= 2:
            variety = parts[0]
            replicate = parts[1]
            if replicate.isdigit():
                return variety, replicate
    except Exception:
        pass
    
    print(f"警告: 无法从文件名 '{filename}' 中按“品种_重复数_...”规则解析，已跳过。")
    return None, None


def parse_bismark_report(filepath):
    """
    解析 Bismark 报告文件，提取所有 key: value 格式的统计信息。
    V3.1更新: 保留原始值，不再去掉 '%' 符号。

    Args:
        filepath (str): Bismark 报告文件的完整路径。

    Returns:
        dict: 包含所有统计项的字典。
    """
    stats = {}
    try:
        real_filepath = os.path.realpath(filepath)
        with open(real_filepath, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split(':\t', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    # 直接获取原始值，不做任何处理
                    value = parts[1].strip()
                    stats[key] = value
    except FileNotFoundError:
        print(f"警告: 找不到文件或链接已损坏 '{filepath}'")
    except Exception as e:
        print(f"警告: 解析文件 '{filepath}' 时发生错误: {e}")
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="统计 WGBS Bismark 比对报告。V3.1: 保留了百分比符号(%)。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("directory", help="必填: 要搜索的 Bismark 报告所在的根目录路径。")
    parser.add_argument("suffix", help="必填: 指定要查找的文件结尾，例如 'bismark_bt2_PE_report.txt'。")
    parser.add_argument("output_file", help="必填: 输出的统计结果表格文件名 (例如 wgbs_summary.tsv)。")
    parser.add_argument("log_file", help="必填: 用于记录所有被处理的文件路径列表的文件名 (例如 processed_files.log)。")
    args = parser.parse_args()

    # 预先定义表头顺序，与Bismark报告完全一致
    DEFINED_HEADER_ORDER = [
        "Sequence pairs analysed in total",
        "Number of paired-end alignments with a unique best hit",
        "Mapping efficiency",
        "Sequence pairs with no alignments under any condition",
        "Sequence pairs did not map uniquely",
        "Sequence pairs which were discarded because genomic sequence could not be extracted",
        "CT/GA/CT",
        "GA/CT/CT",
        "GA/CT/GA",
        "CT/GA/GA",
        "Number of alignments to (merely theoretical) complementary strands being rejected in total",
        "Total number of C's analysed",
        "Total methylated C's in CpG context",
        "Total methylated C's in CHG context",
        "Total methylated C's in CHH context",
        "Total methylated C's in Unknown context",
        "Total unmethylated C's in CpG context",
        "Total unmethylated C's in CHG context",
        "Total unmethylated C's in CHH context",
        "Total unmethylated C's in Unknown context",
        "C methylated in CpG context",
        "C methylated in CHG context",
        "C methylated in CHH context",
        "C methylated in unknown context (CN or CHN)"
    ]

    # --- 1. 查找文件并解析 ---
    all_results = []
    processed_files_log = []

    print(f"正在目录 '{args.directory}' 及其子目录中搜索 *{args.suffix} 文件 (将跟随符号链接)...")

    candidate_files = []
    for root, _, files in os.walk(args.directory, followlinks=True):
        for filename in files:
            if filename.endswith(args.suffix):
                candidate_files.append(os.path.join(root, filename))

    if not candidate_files:
        print("未找到任何匹配后缀的文件。程序退出。")
        return

    for filepath in candidate_files:
        filename = os.path.basename(filepath)
        
        variety, replicate = parse_bismark_filename(filename)
        if not variety:
            continue
            
        stats_data = parse_bismark_report(filepath)
        if not stats_data:
            print(f"警告: 文件 '{filepath}' 为空或无法解析，已跳过。")
            continue
            
        result_entry = {
            'variety': variety,
            'replicate': replicate,
            'stats': stats_data
        }
        all_results.append(result_entry)
        processed_files_log.append(filepath)

    print(f"搜索完成。共找到 {len(candidate_files)} 个匹配文件，其中 {len(all_results)} 个被成功解析。")

    if not all_results:
        print("没有可处理的数据。程序退出。")
        return

    # --- 2. 准备输出 ---
    all_results.sort(key=lambda x: (x['variety'], int(x['replicate'])))
    final_header = ['Variety', 'Replicate'] + DEFINED_HEADER_ORDER

    # --- 3. 写入统计结果文件 ---
    try:
        with open(args.output_file, 'w', encoding='utf-8') as f_out:
            f_out.write('\t'.join(final_header) + '\n')
            
            for result in all_results:
                row_data = [result['variety'], result['replicate']]
                for key in DEFINED_HEADER_ORDER:
                    row_data.append(result['stats'].get(key, 'N/A'))
                
                f_out.write('\t'.join(map(str, row_data)) + '\n')
        
        print(f"\n统计结果已保存到文件: {args.output_file}")

    except IOError as e:
        print(f"错误: 无法写入统计结果文件 '{args.output_file}'。原因: {e}")
        
    # --- 4. 写入日志文件 ---
    try:
        with open(args.log_file, 'w', encoding='utf-8') as f_log:
            f_log.write("# List of files successfully parsed and included in the summary\n")
            f_log.write(f"# Total processed files: {len(processed_files_log)}\n\n")
            processed_files_log.sort()
            for fpath in processed_files_log:
                f_log.write(fpath + '\n')
        print(f"所有被处理文件的路径列表已保存到: {args.log_file}")
    except IOError as e:
        print(f"错误: 无法写入日志文件 '{args.log_file}'。原因: {e}")


if __name__ == "__main__":
    main()
