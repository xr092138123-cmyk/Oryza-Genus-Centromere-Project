#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys

def add_prefix_to_file(filepath):
    """
    将文件名的前缀作为新的一列添加到文件的每一行，并原地保存。

    参数:
        filepath (str): 需要修改的文件的路径。
    """
    # 1. 检查文件是否存在
    if not os.path.exists(filepath):
        print(f"错误: 找不到输入文件 '{filepath}'", file=sys.stderr)
        sys.exit(1)

    # 2. 从文件路径中提取文件名前缀
    # os.path.basename() -> "repeat_region_6617.time"
    # os.path.splitext() -> ("repeat_region_6617", ".time")
    base_name = os.path.basename(filepath)
    prefix, _ = os.path.splitext(base_name)

    print(f"正在处理文件: {filepath}")
    print(f"提取到的前缀是: {prefix}")

    # 3. 将文件所有内容读入内存
    try:
        with open(filepath, 'r') as f:
            original_lines = f.readlines()
    except IOError as e:
        print(f"错误: 无法读取文件 '{filepath}'. 原因: {e}", file=sys.stderr)
        sys.exit(1)
        
    # 4. 构建新的内容
    new_lines = []
    for line in original_lines:
        # 去除行尾的换行符和空白，以防万一
        stripped_line = line.strip()
        if stripped_line:  # 确保不是空行
            # 构建新行： "前缀\t原始行内容\n"
            new_line = f"{prefix}\t{stripped_line}\n"
            new_lines.append(new_line)

    # 5. 将修改后的内容写回原文件
    try:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
    except IOError as e:
        print(f"错误: 无法写入文件 '{filepath}'. 原因: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"成功！文件 '{filepath}' 已被原地修改。")


def main():
    """主执行函数，处理命令行参数"""
    parser = argparse.ArgumentParser(
        description="将文件名的前缀作为第一列添加到文件中，并原地修改文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--input', 
        required=True, 
        help="需要被修改的输入文件路径。"
    )
    args = parser.parse_args()

    add_prefix_to_file(args.input)


if __name__ == "__main__":
    main()
