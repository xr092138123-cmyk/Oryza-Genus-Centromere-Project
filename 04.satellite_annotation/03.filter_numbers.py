import os
import glob

def count_records_in_bed(filepath):
    """
    计算一个BED文件中的记录数量（即行数）。
    """
    count = 0
    try:
        # 使用'with'语句可以确保文件被正确关闭
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            # 遍历文件中的每一行，这比一次性读入内存更高效
            for _ in f:
                count += 1
    except Exception as e:
        print(f"  -> 错误: 无法读取文件 {filepath}: {e}")
        return -1  # 返回一个错误代码，表示读取失败
    return count

def process_and_clean_bed_files():
    """
    主函数：查找、检查并根据记录数清理当前目录下的BED文件。
    """
    # 获取当前工作目录
    current_directory = os.getcwd()
    print(f"开始扫描目录: {current_directory}\n")

    # 使用 glob 查找所有以 .bed 结尾的文件（不区分大小写）
    bed_files = glob.glob('*.bed') + glob.glob('*.BED')
    
    # 去除可能存在的重复项并排序
    bed_files = sorted(list(set(bed_files)))

    if not bed_files:
        print("未在当前目录中找到任何 *.bed 文件。")
        return

    print(f"找到了 {len(bed_files)} 个BED文件，开始检查...\n" + "-"*30)

    # 遍历所有找到的BED文件
    for filename in bed_files:
        filepath = os.path.join(current_directory, filename)

        # 1. 计算文件中的记录数
        record_count = count_records_in_bed(filepath)

        # 如果文件读取失败，则跳过此文件
        if record_count == -1:
            print("-" * 30)
            continue

        print(f"正在检查文件: '{filename}' ... 包含 {record_count} 条记录。")

        # 2. 判断记录数是否小于20，并执行相应操作
        if record_count < 20:
            try:
                os.remove(filepath)
                print(f"  -> 已删除. (原因: 记录数 {record_count} < 20)\n")
            except OSError as e:
                print(f"  -> 删除失败! 错误: {e}\n")
        else:
            print(f"  -> 已保留. (原因: 记录数 {record_count} >= 20)\n")
        
        print("-" * 30)

    print("所有文件检查完毕。")

# 当该脚本被直接运行时，执行主函数
if __name__ == "__main__":
    process_and_clean_bed_files()
