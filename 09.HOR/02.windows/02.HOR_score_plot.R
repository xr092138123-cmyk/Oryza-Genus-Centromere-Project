#!/usr/bin/env Rscript

# --- 加载必要的库 ---
# suppressPackageStartupMessages 用于保持命令行输出整洁
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales)) # 用于坐标轴格式化

# --- 1. 定义命令行参数 ---
parser <- ArgumentParser(description = "从TSV文件绘制HOR score的基因组分布柱状图。")
parser$add_argument("--input", required = TRUE, help = "输入的TSV文件路径。")
parser$add_argument("--output", required = TRUE, help = "输出的图片文件路径 (例如: hor_plot.png, hor_plot.pdf)。")

args <- parser$parse_args()

# --- 2. 读取数据与错误检查 ---
if (!file.exists(args$input)) {
  stop(paste("错误：找不到输入文件", args$input))
}

cat("正在读取文件:", args$input, "...\n")
# 读取TSV，不自动转换字符为因子
df <- read_tsv(args$input, col_types = cols(), show_col_types = FALSE)

# 检查必需列
required_cols <- c("label", "HOR_score_percent")
if (!all(required_cols %in% colnames(df))) {
  stop(paste("错误：输入文件必须包含以下列:", paste(required_cols, collapse = ", ")))
}

if (nrow(df) == 0) {
  stop("错误：输入文件为空。")
}

# --- 3. 数据处理 ---
cat("正在处理数据...\n")

plot_data <- df %>%
  # 1. 解析 'label' 列 (格式如 Chr:start-end)
  # 使用正则表达式提取 start 和 end。
  # convert = TRUE 会自动将提取出的数字字符串转换为数值型
  tidyr::extract(
    col = label,
    into = c("chrom", "start", "end"),
    regex = "^(.*?):([0-9]+)-([0-9]+)$",
    convert = TRUE,
    remove = FALSE # 保留原label列以防万一
  ) %>%
  # 2. 过滤掉解析失败的行 (如果有)
  filter(!is.na(start) & !is.na(end)) %>%
  # 3. 确保 HOR_score_percent 是数值型
  mutate(HOR_score_percent = as.numeric(HOR_score_percent)) %>%
  # 4. 筛选分数大于0的数据 (对应 Python 脚本中的逻辑)
  filter(HOR_score_percent > 0)

# 检查处理后是否还有数据
if (nrow(plot_data) == 0) {
  cat("提示：没有找到任何 HOR score > 0 且格式正确的记录，无法生成图形。\n")
  quit(status = 0)
}

cat(sprintf("成功处理了 %d 条记录。\n", nrow(plot_data)))

# --- 4. 绘图 ---
cat("正在生成图形...\n")

# 使用 geom_rect 替代 barplot
# Python的 ax.bar(align='edge', width=width) 等同于画一个从 start 到 end 的矩形
p <- ggplot(plot_data) +
  geom_rect(aes(
    xmin = start, 
    xmax = end, 
    ymin = 0, 
    ymax = HOR_score_percent
  ), fill = "#1E90FF", alpha = 0.9) +
  
  # 设置坐标轴标签和标题
  labs(
    title = "Distribution of HOR Scores along the Chromosome",
    x = "Genomic Position",
    y = "HOR Score (%)"
  ) +
  
  # 格式化 X 轴为 Mb (Megabases)
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = " Mb")) +
  
  # 设置主题
  theme_classic() +
  

  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 10, color = "black")
  )

# --- 5. 保存图形 ---
cat("正在保存图形到:", args$output, "...\n")

tryCatch({
  ggsave(
    filename = args$output, 
    plot = p, 
    width = 16, 
    height = 6, 
    dpi = 300,
    limitsize = FALSE
  )
  cat("图形已成功保存。\n")
}, error = function(e) {
  stop(paste("无法保存图形文件:", e$message))
})
