#!/usr/bin/env python3

import argparse
import sys

# --- 配置分类规则 ---
# 将规则定义在这里，方便未来修改或添加新类型
# 格式为: {"分类名称": 典型长度}
CLASSIFICATION_RULES = {
    "CEN103": 103,
    "CEN126": 126,
    "CEN155": 155,
    "CEN165": 165,
    "CEN180": 180,
}
# 定义长度容差
TOLERANCE = 5

def get_classification(length: int) -> str:
    """
    根据给定的长度和预设规则返回分类名称。

    Args:
        length (int): 记录的长度。

    Returns:
        str: 分类名称 (e.g., 'CEN155') 或 'atypical'。
    """
    # 遍历所有规则
    for name, typical_len in CLASSIFICATION_RULES.items():
        # 检查长度是否在 [typical - 5, typical + 5) 范围内
        if (typical_len - TOLERANCE) <= length < (typical_len + TOLERANCE):
            return name  # 找到匹配项，立即返回分类名称
    
    # 如果循环结束都没有找到匹配项，则返回 'atypical'
    return "atypical"

def process_file(input_file: str, output_file: str):
    """
    处理输入文件，进行分类，并写入输出文件。

    Args:
        input_file (str): 输入的TSV文件路径。
        output_file (str): 输出的TSV文件路径。
    """
    print(f"开始处理文件: {input_file}")
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # 逐行读取输入文件
            for line_num, line in enumerate(infile, 1):
                # 跳过空行
                if not line.strip():
                    continue
                
                # TSV文件是tab分隔的
                parts = line.strip().split('\t')
                
                # --- 健壮性检查 ---
                # 检查行是否有足够的列
                if len(parts) < 7:
                    print(f"警告 (行 {line_num}): 列数少于7，将原样写入。内容: '{line.strip()}'")
                    outfile.write(line)
                    continue
                
                try:
                    # 提取第七列的长度并转换为整数
                    record_length = int(parts[6])
                except ValueError:
                    print(f"警告 (行 {line_num}): 第七列不是有效整数，将原样写入。内容: '{line.strip()}'")
                    outfile.write(line)
                    continue

                # --- 核心逻辑 ---
                # 1. 获取分类名称
                new_type = get_classification(record_length)
                
                # 2. 替换第四列的内容 (parts[3]因为列表索引从0开始)
                parts[3] = new_type
                
                # 3. 将修改后的列表重新组合成一个tab分隔的字符串
                output_line = "\t".join(parts) + "\n"
                
                # 4. 写入输出文件
                outfile.write(output_line)
                
        print(f"处理完成！结果已成功保存到: {output_file}")

    except FileNotFoundError:
        print(f"错误：找不到输入文件 '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理文件时发生未知错误: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """
    主函数，用于设置和解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="根据第七列的长度对TSV文件进行分类，并更新第四列的标识符。"
    )
    
    parser.add_argument('--input', '-i',
                        required=True,
                        help="输入的TSV文件路径。")
                        
    parser.add_argument('--output', '-o',
                        required=True,
                        help="输出的分类后的TSV文件路径。")
                        
    args = parser.parse_args()
    
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
