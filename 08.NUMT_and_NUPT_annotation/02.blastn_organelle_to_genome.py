import argparse
import subprocess
import os
import sys

def run_blast_pipeline(query_file, target_file):
    """
    执行 makeblastdb 和 blastn 的完整流程，并将数据库文件整理到单独的文件夹中。

    Args:
        query_file (str): 查询FASTA文件的路径。
        target_file (str): 目标FASTA文件的路径。
    """
    print("--- 开始执行BLAST流程 ---")

    # --- 步骤 0: 检查输入文件和依赖项是否存在 ---
    if not os.path.exists(query_file):
        print(f"错误: 查询文件 '{query_file}' 未找到！", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(target_file):
        print(f"错误: 目标文件 '{target_file}' 未找到！", file=sys.stderr)
        sys.exit(1)

    # --- 步骤 1: 准备文件名和路径---
    # 从目标文件名自动生成文件夹名 (例如 NIP-T2T.fa -> NIP-T2T)
    target_basename = os.path.basename(target_file)
    db_folder = os.path.splitext(target_basename)[0]

    # 构建数据库文件的完整路径前缀。数据库文件将位于 db_folder 中，且文件名与文件夹名相同。
    # 例如: NIP-T2T/NIP-T2T
    db_path_prefix = os.path.join(db_folder, db_folder)
    
    # 输出文件名逻辑保持不变，结果仍在当前目录
    query_base = os.path.splitext(os.path.basename(query_file))[0]
    output_file = f"{query_base}.{db_folder}.blast.txt"

    print(f"查询文件: {query_file}")
    print(f"目标文件: {target_file}")
    print(f"数据库文件夹将是: {db_folder}/")
    print(f"数据库路径前缀将是: {db_path_prefix}")
    print(f"输出文件将是: {output_file}")
    print("-" * 20)

    # --- 步骤 2: 创建数据库文件夹 ---
    print(f"\n[1/3] 正在创建数据库文件夹...")
    try:
        os.makedirs(db_folder, exist_ok=True)
        print(f"文件夹 '{db_folder}' 已创建或已存在。")
    except OSError as e:
        print(f"错误: 创建文件夹 '{db_folder}' 失败: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 步骤 3: 执行 makeblastdb 命令---
    print("\n[2/3] 正在构建BLAST数据库...")
    makeblastdb_cmd = [
        'makeblastdb',
        '-in', target_file,
        '-dbtype', 'nucl',
        '-parse_seqids',
        '-out', db_path_prefix  # <--- 使用包含文件夹的路径
    ]
    print(f"执行命令: {' '.join(makeblastdb_cmd)}")
    
    try:
        subprocess.run(makeblastdb_cmd, check=True, capture_output=True, text=True)
        print("数据库构建成功！")
    except FileNotFoundError:
        print("\n错误: 'makeblastdb' 命令未找到。", file=sys.stderr)
        print("请确保 NCBI BLAST+ 已正确安装，并且其路径已添加到您的系统环境变量 (PATH) 中。", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\n错误: 'makeblastdb' 执行失败。返回码: {e.returncode}", file=sys.stderr)
        print(f"错误信息:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

    # --- 步骤 4: 执行 blastn 命令 ---
    print("\n[3/3] 正在运行 blastn 比对...")
    blastn_cmd = [
        'blastn',
        '-query', query_file,
        '-db', db_path_prefix,  
        '-out', output_file,
        '-outfmt', '6'
    ]
    print(f"执行命令: {' '.join(blastn_cmd)}")

    try:
        subprocess.run(blastn_cmd, check=True)
        print("blastn 比对成功！")
    except FileNotFoundError:
        print("\n错误: 'blastn' 命令未找到。", file=sys.stderr)
        print("请确保 NCBI BLAST+ 已正确安装，并且其路径已添加到您的系统环境变量 (PATH) 中。", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\n错误: 'blastn' 执行失败。返回码: {e.returncode}", file=sys.stderr)
        print(f"错误信息:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)
    
    print("\n--- 流程执行完毕 ---")
    print(f"比对结果已保存到: {output_file}")


def main():
    """
    主函数，用于解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description="一个自动化运行 makeblastdb 和 blastn 的Python脚本，并将数据库文件整理到文件夹中。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-q", "--query", 
        required=True, 
        help="[必需] 查询序列文件 (FASTA格式)。"
    )
    parser.add_argument(
        "-t", "--target", 
        required=True, 
        help="[必需] 目标/参考序列文件 (FASTA格式)，用于构建数据库。"
    )
    
    args = parser.parse_args()
    
    run_blast_pipeline(args.query, args.target)

if __name__ == "__main__":
    main()
