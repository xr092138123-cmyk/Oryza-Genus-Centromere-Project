#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

def load_time_data(time_filepath):
    """
    读取时间文件，并将其内容加载到一个字典中以便快速查找。
    键是第一列的ID，值是最后一列的时间。

    参数:
        time_filepath (str): time文件的路径。

    返回:
        dict: 一个字典，键是记录ID (第一列)，值是插入时间 (最后一列)。
    """
    time_map = {}
    print(f"正在加载时间文件: {time_filepath}...")
    try:
        with open(time_filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # 使用split()处理可能存在的多个空格或制表符
                columns = line.split()
                if len(columns) >= 2: # 至少需要一个ID和时间
                    record_id = columns[0]
                    insertion_time = columns[-1] # -1 表示最后一列
                    time_map[record_id] = insertion_time
                else:
                    print(f"警告: 跳过时间文件中格式不正确的行: '{line}'", file=sys.stderr)
    except FileNotFoundError:
        print(f"错误: 找不到时间文件 '{time_filepath}'", file=sys.stderr)
        sys.exit(1)
        
    print(f"成功从时间文件中加载 {len(time_map)} 条记录。")
    return time_map

def merge_files(bed_filepath, time_map, output_filepath):
    """
    逐行读取BED文件，从第五列匹配ID，合并时间数据，并写入输出文件。

    参数:
        bed_filepath (str): bed文件的路径。
        time_map (dict): 包含时间和ID映射的字典。
        output_filepath (str): 输出文件的路径。
    """
    not_found_count = 0
    total_lines = 0
    print(f"正在处理BED文件: {bed_filepath}...")
    
    try:
        with open(bed_filepath, 'r') as bed_file, open(output_filepath, 'w') as out_file:
            for line in bed_file:
                original_line = line.strip()
                if not original_line:
                    continue
                
                total_lines += 1
                bed_columns = original_line.split('\t')
                
                # 关键点：ID在第五列 (索引为4)
                if len(bed_columns) >= 5:
                    bed_id = bed_columns[4]
                    
                    # 使用.get()方法查找时间，如果找不到，则返回默认值 'NA'
                    insertion_time = time_map.get(bed_id, 'NA')
                    
                    if insertion_time == 'NA':
                        not_found_count += 1
                    
                    # 构建新的输出行
                    new_line = f"{original_line}\t{insertion_time}\n"
                    out_file.write(new_line)
                else:
                    # 如果行格式不正确，直接写入原始行并附上NA
                    out_file.write(f"{original_line}\tNA\n")
                    print(f"警告: BED文件中发现格式不正确的行: '{original_line}'", file=sys.stderr)

    except FileNotFoundError:
        print(f"错误: 找不到BED文件 '{bed_filepath}'", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"错误: 无法写入输出文件 '{output_filepath}'. 原因: {e}", file=sys.stderr)
        sys.exit(1)
        
    print("\n处理完成！")
    print(f"总共处理了 {total_lines} 行BED记录。")
    if not_found_count > 0:
        print(f"警告: 有 {not_found_count} 条记录在时间文件中未找到匹配项，已在输出中标注为 'NA'。")


def main():
    """主执行函数"""
    parser = argparse.ArgumentParser(
        description="将TE插入时间数据合并到BED文件的末尾。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--bed', 
        required=True, 
        help="输入的TE BED文件路径 (ID在第五列)。"
    )
    parser.add_argument(
        '--time', 
        required=True, 
        help="包含插入时间信息的TSV文件路径 (ID在第一列，时间在最后一列)。"
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help="输出的合并后BED文件的路径。"
    )
    args = parser.parse_args()

    # 1. 加载时间数据到字典
    time_data = load_time_data(args.time)

    # 2. 合并文件
    merge_files(args.bed, time_data, args.output)
    
    print(f"\n成功！合并后的文件已保存至: '{args.output}'")


if __name__ == "__main__":
    main()
