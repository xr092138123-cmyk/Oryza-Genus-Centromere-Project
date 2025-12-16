import pandas as pd
import argparse
import sys
import subprocess
import os

def check_bedtools_installed():
    """检查 bedtools 是否已安装并且可用。"""
    try:
        # 运行 'bedtools --version'，并将输出重定向到DEVNULL以保持终端干净
        subprocess.run(['bedtools', '--version'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def run_bedtools_summary(input_file, output_summary_file):
    """
    使用 bedtools merge 合并区间，然后用 pandas 统计每个染色体的总长度。

    Args:
        input_file (str): 输入的BED文件名。
        output_summary_file (str): 输出的染色体长度统计文件名。
    """
    # --- 步骤 0: 检查依赖 ---
    if not check_bedtools_installed():
        print("错误: 'bedtools' 未找到或无法执行。", file=sys.stderr)
        print("请确保 bedtools 已正确安装，并且其路径已添加到您的系统环境变量 (PATH) 中。", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(input_file):
        print(f"错误: 输入文件 '{input_file}' 未找到。", file=sys.stderr)
        sys.exit(1)
        
    # 定义一个临时的中间文件名
    merged_bed_file = f"{os.path.splitext(input_file)[0]}.merged.tmp.bed"

    try:
        # --- 步骤 1: 调用 bedtools merge ---
        print(f"[1/4] 正在使用 bedtools merge 合并区间...")
        bedtools_cmd = ['bedtools', 'merge', '-i', input_file]
        
        # 打开一个文件句柄用于写入 bedtools 的输出
        with open(merged_bed_file, 'w') as f_out:
            # 执行命令，将标准输出重定向到我们打开的文件
            subprocess.run(bedtools_cmd, stdout=f_out, check=True)
        print(f"-> 合并后的区间已临时保存到: {merged_bed_file}")

        # --- 步骤 2: 使用 pandas 读取和处理合并后的文件 ---
        print(f"\n[2/4] 正在读取合并后的文件并计算长度...")
        # 定义合并后文件的列名
        merged_columns = ['chrom', 'start', 'end']
        df_merged = pd.read_csv(
            merged_bed_file, 
            sep='\t', 
            header=None, 
            names=merged_columns
        )
        
        # 计算每个合并后区间的长度
        df_merged['length'] = df_merged['end'] - df_merged['start']
        
        # 按染色体分组并对长度求和
        length_stats = df_merged.groupby('chrom')['length'].sum().reset_index()
        
        # 重命名列以符合最终输出要求
        length_stats.columns = ['chromosome', 'total_length']
        print("-> 长度计算完成。")

        # --- 步骤 3: 保存最终的统计结果 ---
        print(f"\n[3/4] 正在将统计结果保存到: {output_summary_file} ...")
        length_stats.to_csv(
            output_summary_file, 
            sep='\t', 
            header=True, # 包含标题行，更具可读性
            index=False
        )
        print("\n处理成功！")
        print(f"统计文件已保存至: {output_summary_file}")

    except subprocess.CalledProcessError as e:
        print(f"错误: 'bedtools merge' 执行失败。返回码: {e.returncode}", file=sys.stderr)
        print(f"请检查输入文件 '{input_file}' 是否为有效的BED格式。", file=sys.stderr)
    except Exception as e:
        print(f"处理过程中发生错误: {e}", file=sys.stderr)
    finally:
        # --- 步骤 4: 清理临时文件 ---
        if os.path.exists(merged_bed_file):
            print(f"\n[4/4] 正在清理临时文件: {merged_bed_file} ...")
            os.remove(merged_bed_file)


def main():
    """
    主函数，用于解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="使用 bedtools merge 合并BED区间并统计每个染色体的总覆盖长度。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--input_file", 
        help="输入的BED文件名。\n例如: my_regions.bed"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="输出的长度统计文件名 (两列，制表符分隔)。\n例如: coverage_summary.txt"
    )
    
    args = parser.parse_args()
    run_bedtools_summary(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
