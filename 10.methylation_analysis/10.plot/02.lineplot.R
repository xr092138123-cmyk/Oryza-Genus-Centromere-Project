# 加载包（确保已安装 tidyverse、fuzzyjoin）
library(tidyverse)
library(fuzzyjoin)  


# 1. 读取样本-稻种映射（空格分隔）
name_mapping <- read.table(
  "name.txt", 
  sep = "",  # 自动识别空格/制表符
  header = FALSE, 
  col.names = c("sample", "species"),
  stringsAsFactors = FALSE
)

target_sample <- "HZ"
target_species <- name_mapping %>% 
  filter(sample == target_sample) %>% 
  pull(species)


# 2. 读取着丝粒区域（制表符分隔）
centromere <- read.delim(
  "Oryza13_centromere_pos.xls", 
  sep = "\t", 
  header = TRUE,
  stringsAsFactors = FALSE
)

# 筛选目标稻种+hap1，并重命名染色体列避免冲突
centro_hap1 <- centromere %>% 
  filter(Name == target_species, Hap == "hap1") %>% 
  select(Chr, Start, End) %>% 
  rename(centro_chr = Chr)  # 重命名为centro_chr


# 3. 读取甲基化数据（2kb区间）
methylation <- read.csv(
  "HZ.hap1.5x.csv", 
  stringsAsFactors = FALSE
)

# 过滤无效数据 + 确保chr为字符型
methyl_valid <- methylation %>% 
  filter(total_coverage > 0) %>% 
  drop_na(avg_methylation) %>% 
  mutate(
    across(c(start, end), as.integer),
    chr = as.character(chr)  # 避免因子类型冲突
  )


# 4. 匹配“甲基化区间”与“着丝粒区域”（区间重叠判断）
methyl_centro <- methyl_valid %>% 
  fuzzy_left_join(
    centro_hap1, 
    by = c(
      "chr" = "centro_chr",  # 染色体匹配
      "start" = "End",       # 甲基化start < 着丝粒End
      "end" = "Start"        # 甲基化end > 着丝粒Start
    ),
    match_fun = list(`==`, `<`, `>`)
  ) %>% 
  filter(!is.na(Start)) %>%  # 仅保留着丝粒内的区间
  select(-c(Start, End))     # 移除冗余字段


# 5. 整理数据：按染色体+区间起始排序
methyl_centro_sorted <- methyl_centro %>% 
  arrange(chr, start) %>%  # 先按染色体，再按区间起始
  group_by(chr) %>%        # 按染色体分组计算相对位置
  mutate(
    rel_pos = start - min(start),        # 染色体内部相对位置（bp）
    rel_pos_kb = rel_pos / 1000          # 转换为kb
  ) %>% 
  ungroup()


# 6. 绘制优化版折线图
plot <- ggplot(
  methyl_centro_sorted, 
  aes(x = rel_pos_kb, y = avg_methylation, color = context)
) +
  # 细化线条（linewidth=0.5）+ 缩小点（size=2）
  geom_line(linewidth = 0.2) +        
  geom_point(size = 0.2, alpha = 0.7) +  # 半透明避免重叠
  facet_wrap(
    vars(chr),        # 每个染色体一个子图
    scales = "free_x",# x轴自由缩放（适配不同染色体长度）
    ncol = 1          # 垂直排列子图
  ) +
  labs(
    x = "Relative Position in Centromere (kb)",
    y = "Methylation Level (%)",
    title = paste0(target_sample, " (", target_species, ")"),
    color = "Context"
  ) +
  # 让y轴更舒展（保留0-100范围，可根据数据调整）
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # x轴文字缩小
    axis.text.y = element_text(size = 10),                        # 增大y轴文字
    strip.text = element_text(size = 11, face = "bold"),          # 子图标题加粗
    legend.position = "right",                                    # 图例放右侧
    legend.text = element_text(size = 10)                         # 图例文字大小
  )



ggsave(
  "centromere_methylation_optimized.5x.png", 
  plot, 
  width = 10, 
  height = 40,  # 垂直排列时增大高度
  dpi = 300
)
