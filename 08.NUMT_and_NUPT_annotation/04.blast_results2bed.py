import argparse
import sys

def blast_to_bed(blast_file, bed_file):
    """
    将 BLAST outfmt 6 结果文件转换为 BED 格式文件。

    Args:
        blast_file (str): 输入的 BLAST (outfmt 6) 文件名。
        bed_file (str): 输出的 BED 文件名。
    """
    print(f"正在读取 BLAST 文件: {blast_file}")
    
    lines_processed = 0
    lines_skipped = 0
    
    try:
        with open(blast_file, 'r') as infile, open(bed_file, 'w') as outfile:
            for line in infile:
                # 跳过空行或注释行
                if not line.strip() or line.startswith('#'):
                    continue
                
                try:
                    # 分割行，outfmt 6 是制表符或空格分隔的
                    fields = line.strip().split()
                    
                    # 确保至少有12列
                    if len(fields) < 12:
                        lines_skipped += 1
                        continue

                    # --- 1. 从BLAST结果中提取所需信息 ---
                    # qseqid (查询序列ID) -> BED name
                    qseqid = fields[0]
                    # sseqid (目标序列ID) -> BED chrom
                    sseqid = fields[1]
                    # sstart, send (目标序列上的比对起止位置, 1-based)
                    sstart = int(fields[8])
                    send = int(fields[9])
                    # bitscore -> BED score
                    bitscore = fields[11]

                    # --- 2. 转换为BED格式的字段 ---
                    # BED chrom: 就是目标序列的ID
                    bed_chrom = sseqid
                    
                    # BED name: 使用查询序列的ID
                    bed_name = qseqid
                    
                    # BED score: 使用bitscore
                    bed_score = bitscore
                    
                    # BED strand:
                    # 如果 sstart < send, 说明是正向链 '+'
                    # 如果 sstart > send, 说明是反向链 '-'
                    if sstart < send:
                        bed_strand = '+'
                        # BED chromStart 是 0-based，所以要减1
                        bed_start = sstart - 1
                        bed_end = send
                    else: # sstart > send
                        bed_strand = '-'
                        # 即使在反向链，BED的start也总是小于end
                        bed_start = send - 1
                        bed_end = sstart
                    
                    # --- 3. 格式化并写入BED文件 ---
                    # BED格式是 tab 分隔的
                    bed_line = f"{bed_chrom}\t{bed_start}\t{bed_end}\t{bed_name}\t{bed_score}\t{bed_strand}\n"
                    outfile.write(bed_line)
                    
                    lines_processed += 1
                    
                except (ValueError, IndexError) as e:
                    # 如果某一行格式不正确（例如，坐标不是数字），则跳过并打印警告
                    print(f"警告: 跳过格式错误的行: '{line.strip()}' -> 错误: {e}", file=sys.stderr)
                    lines_skipped += 1
                    continue
                    
    except FileNotFoundError:
        print(f"错误: 输入文件 '{blast_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
        
    print("\n转换完成！")
    print(f"成功处理了 {lines_processed} 行。")
    if lines_skipped > 0:
        print(f"跳过了 {lines_skipped} 行格式错误的行。")
    print(f"BED 文件已保存到: {bed_file}")

def main():
    """
    主函数，用于解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="将 BLASTn outfmt 6 结果文件转换为标准的6列BED文件。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--blast_file", 
        help="输入的BLAST结果文件 (outfmt 6 格式)。\n例如: blast.txt"
    )
    parser.add_argument(
        '-o', "--bed_file", 
        help="输出的BED文件名。\n例如: output.bed"
    )
    
    args = parser.parse_args()
    
    blast_to_bed(args.blast_file, args.bed_file)

if __name__ == "__main__":
    main()
