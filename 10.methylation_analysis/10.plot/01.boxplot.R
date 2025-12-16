library(tidyverse)

# ----------------------
# 关键修复：设置兼容所有系统的字体
# ----------------------
theme_set(theme_bw() + theme(text = element_text(family = "sans")))

# ----------------------
# 1. 路径设置
# ----------------------
name_file <- "name.txt"  # 确保该文件在脚本运行目录，或填写绝对路径
centromere_file <- "/share/org/YZWL/yzwl_liuyy/wangdan/WGBS/Oryza13_centromere_pos.xls"
base_dir <- "/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/Result/result_hap1/2kb"  # 所有主文件夹的上级目录
output_prefix <- "methylation_"

# ----------------------
# 2. 读取name.txt
# ----------------------
name_mapping <- read.table(
  name_file,
  sep = "",
  header = FALSE,
  col.names = c("file_prefix", "centro_name"),
  fill = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE
)

if (ncol(name_mapping) != 2) {
  stop("name.txt parsing error! Ensure it has two columns separated by spaces")
}

# ----------------------
# 3. 读取着丝粒文件
# ----------------------
centromere_df <- read_delim(
  centromere_file,
  delim = "\t",
  col_names = TRUE,
  trim_ws = TRUE,
  show_col_types = FALSE
) %>% 
  rename(Name = Name, Hap = Hap, Chr = Chr, Start = Start, End = End) %>% 
  mutate(Start = as.integer(Start), End = as.integer(End))

required_centro_cols <- c("Name", "Hap", "Chr", "Start", "End")
if (!all(required_centro_cols %in% colnames(centromere_df))) {
  missing_cols <- setdiff(required_centro_cols, colnames(centromere_df))
  stop(paste("Centromere file missing columns:", paste(missing_cols, collapse = ", ")))
}

# ----------------------
# 4. 处理甲基化文件
# ----------------------
# 递归搜索所有子文件夹中符合格式的CSV文件（如HZ.hap1.1.csv）
methyl_files <- list.files(
  path = base_dir,
  pattern = "^.*\\.hap\\d+\\.\\d+\\.csv$",  # 匹配改名后的文件名格式
  full.names = TRUE,
  recursive = TRUE  # 关键：递归搜索子文件夹
)

# 初始化列表存储数据
all_data_list <- list()

for (file in methyl_files) {
  file_name <- basename(file)
  
  # 提取物种前缀、单倍型编号（如HZ.hap1.1.csv -> 前缀HZ，hap1，重复1）
  if (str_detect(file_name, "^(.*?)\\.hap(\\d+)\\.(\\d+)\\.csv$")) {
    species_prefix <- str_extract(file_name, "^(.*?)\\.hap(\\d+)\\.(\\d+)\\.csv$", group = 1)
    hap_num <- str_extract(file_name, "^(.*?)\\.hap(\\d+)\\.(\\d+)\\.csv$", group = 2)
    rep_num <- str_extract(file_name, "^(.*?)\\.hap(\\d+)\\.(\\d+)\\.csv$", group = 3)  # 提取重复编号（1/2/3）
    hap <- paste0("hap", hap_num)
  } else {
    warning(paste("文件", file_name, "格式不符，跳过"))
    next
  }
  
  # 匹配name.txt中的物种名称
  species_match <- filter(name_mapping, file_prefix == species_prefix)
  if (nrow(species_match) == 0) {
    warning(paste("前缀", species_prefix, "未在name.txt中找到，跳过文件", file_name))
    next
  }
  species <- species_match$centro_name
  
  # 读取甲基化文件
  methyl_df <- tryCatch({
    read_csv(file, show_col_types = FALSE) %>%
      rename(chr = chr, start = start, end = end, context = context, avg_methylation = avg_methylation) %>%
      filter(!is.na(avg_methylation))
  }, error = function(e) {
    warning(paste("读取文件", file_name, "失败：", e$message))
    return(NULL)
  })
  
  if (is.null(methyl_df) || nrow(methyl_df) == 0) next
  
  # 匹配对应的着丝粒区域
  centro_regions <- filter(centromere_df, Name == species, Hap == hap)
  
  # 标记着丝粒/非着丝粒区域
  methyl_df <- methyl_df %>%
    rowwise() %>%
    mutate(
      in_centro = if(nrow(centro_regions) > 0) {
        any(centro_regions$Chr == chr & !(end < centro_regions$Start | start > centro_regions$End))
      } else {
        FALSE
      },
      region_type = ifelse(in_centro, "Centromere", "Non-Centromere"),
      species = species,
      hap = hap,
      replicate = rep_num  # 新增：记录重复编号
    ) %>%
    ungroup() %>%
    select(species, hap, replicate, context, region_type, methylation = avg_methylation)
  
  all_data_list[[length(all_data_list) + 1]] <- methyl_df
}

# 合并所有数据
all_data <- bind_rows(all_data_list)

if (nrow(all_data) == 0) {
  stop("无有效数据用于绘图，请检查文件路径和格式")
}

# ----------------------
# 5. 绘图（按context生成，图片名包含所有样本信息）
# ----------------------
context_types <- unique(all_data$context)

for (context in context_types) {
  plot_data <- filter(all_data, context == !!context)
  if (nrow(plot_data) == 0) next
  
  # 组合物种、单倍型、重复编号作为x轴标签（如HZ_hap1_1）
  plot_data <- plot_data %>% 
    mutate(sample_id = paste0(species, "_", hap, "_", replicate))
  
  # 绘制箱线图
  plot <- ggplot(plot_data, aes(x = sample_id, y = methylation, fill = region_type)) +
    geom_boxplot(outlier.size = 0.8, width = 0.5, outlier.shape = NA) +
    labs(
      title = paste(context, "Methylation Level"),
      x = "Region Type", 
      y = "Methylation (%)", 
      fill = "Region Type"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      text = element_text(family = "sans")
    )
  
  # 输出图片名：methylation_HZ.hap1.3.CG.pdf 等（按context汇总所有样本）
  output_file <- paste0(output_prefix, context, ".pdf")
  ggsave(
    output_file, 
    plot, 
    width = 10 + n_distinct(plot_data$sample_id) * 0.8,  # 根据样本数量调整宽度
    height = 8, 
    device = "pdf"
  )
  
  cat("已生成", context, "甲基化箱线图：", output_file, "\n")
}

cat("所有图片生成完成！\n")
