
folders=(YZ002.2)
# Bismark路径和参数
BISMARK="bismark_genome_preparation"
BOWTIE2_PATH="bowtie2-2.5.2-linux-x86_64/"

# 循环处理每个基因组文件夹
for folder in "${folders[@]}"; do
    echo -e "\n===== 正在处理基因组: $folder ====="
    
    # 检查文件夹是否存在
    if [ ! -d "$folder" ]; then
        echo "错误: 文件夹 '$folder' 不存在，跳过..."
        continue
    fi
    
    # 执行Bismark索引构建
    $BISMARK --bowtie2 --path_to_aligner "$BOWTIE2_PATH" --verbose "./$folder"
    
    # 检查命令执行状态
    if [ $? -eq 0 ]; then
        echo "正确 基因组 '$folder' 索引构建完成"
    else
        echo "错误 基因组 '$folder' 索引构建失败"
    fi
done

echo -e "\n===== 所有基因组索引构建完成 ====="