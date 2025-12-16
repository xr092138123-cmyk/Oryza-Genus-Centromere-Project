#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from collections import defaultdict

def parse_attributes(attr_string):
    """
    解析BED/GFF文件第七列的键值对属性字符串。
    例如: "id=lLTR_1;parent=repeat_region_1;..."
    
    返回:
        一个属性字典。
    """
    attributes = {}
    parts = attr_string.strip().split(';')
    for part in parts:
        if '=' in part:
            # 只在第一个'='处分割，以防值中也包含'='
            key, value = part.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes

def find_and_write_ltr_pairs(input_filepath, output_dir):
    """
    主处理函数：查找配对的LTR并写入单独的文件。
    """
    # 使用嵌套的defaultdict可以极大地简化代码
    # 结构: {parent_id: {ltr_index: {'lLTR': line, 'rLTR': line}}}
    ltr_data = defaultdict(lambda: defaultdict(dict))
    
    print(f"正在读取和解析输入文件: {input_filepath}...")
    try:
        with open(input_filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                columns = line.split('\t')
                # LTR记录的第五列通常是 'long_terminal_repeat'
                if len(columns) >= 5 and columns[4] == 'long_terminal_repeat':
                    if len(columns) < 7:
                        print(f"警告: 跳过格式不正确的行: {line}", file=sys.stderr)
                        continue
                        
                    attrs = parse_attributes(columns[6])
                    
                    if 'id' in attrs and 'parent' in attrs:
                        ltr_id = attrs['id']
                        parent_id = attrs['parent']
                        
                        # 使用正则表达式匹配 lLTR_xxx 或 rLTR_xxx 格式
                        match = re.match(r'([lr]LTR)_(\S+)', ltr_id)
                        if match:
                            ltr_type = match.group(1) # 'lLTR' or 'rLTR'
                            ltr_index = match.group(2) # '1', '2', etc.
                            
                            # 存储信息
                            ltr_data[parent_id][ltr_index][ltr_type] = line
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{input_filepath}'", file=sys.stderr)
        sys.exit(1)

    print("文件解析完成，开始查找配对并写入文件...")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    pair_count = 0
    # 遍历收集到的数据，查找完整的配对
    for parent_id, ltr_indexes in ltr_data.items():
        for ltr_index, ltr_pair in ltr_indexes.items():
            # 检查是否同时存在 'lLTR' 和 'rLTR'
            if 'lLTR' in ltr_pair and 'rLTR' in ltr_pair:
                pair_count += 1
                
                # 构建输出文件名
                output_filename = os.path.join(output_dir, f"{parent_id}.bed")
                
                try:
                    with open(output_filename, 'w') as out_file:
                        # 写入配对的两条记录
                        out_file.write(ltr_pair['lLTR'] + '\n')
                        out_file.write(ltr_pair['rLTR'] + '\n')
                except IOError as e:
                    print(f"错误: 无法写入文件 '{output_filename}'. 原因: {e}", file=sys.stderr)

    print("\n处理完成！")
    print(f"总共找到了 {pair_count} 对配对的LTR。")
    print(f"输出文件已保存至目录: '{output_dir}'")


def main():
    """主执行函数，处理命令行参数"""
    parser = argparse.ArgumentParser(
        description="从BED文件中查找配对的左右LTR (long_terminal_repeat)，并为每对创建一个新的BED文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help="包含LTR信息的输入BED文件路径。"
    )
    parser.add_argument(
        '--output_dir', 
        required=True, 
        help="用于存放输出的配对LTR文件的目录路径。\n如果目录不存在，脚本将自动创建它。"
    )
    args = parser.parse_args()

    find_and_write_ltr_pairs(args.input, args.output_dir)


if __name__ == "__main__":
    main()
