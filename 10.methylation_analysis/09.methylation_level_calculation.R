library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(ggplot2)

# ======================
# 核心函数：计算2kb窗口甲基化（含数据清洗）
# ======================
calculate_methylation_windows <- function(file_path, context, fai_path = NULL) {
  # 1. 读取数据并严格检查格式
  if (!file.exists(file_path)) stop("文件不存在: ", file_path)
  cat("读取文件:", file_path, "\n")
  
  # 读取数据（强制指定列类型，避免自动转换错误）
  meth_data <- tryCatch({
    fread(
      file_path,
      skip = "chrbase",  # 跳过表头行
      col.names = c("chrbase", "chr", "base", "strand", "coverage", "freqC", "freqT"),
      colClasses = list(
        character = c("chrbase", "chr", "strand"),  # 字符型列
        integer = c("base", "coverage"),            # 整数型列
        numeric = c("freqC", "freqT")               # 数值型列
      ),
      data.table = FALSE  # 返回data.frame
    )
  }, error = function(e) {
    cat("警告：fread读取失败，改用read.table重试\n")
    read.table(
      file_path,
      header = FALSE,
      skip = 1,  # 跳过表头行
      col.names = c("chrbase", "chr", "base", "strand", "coverage", "freqC", "freqT"),
      stringsAsFactors = FALSE
    )
  })
  
  # 2. 严格数据清洗（核心步骤）
  # 检查列数是否正确
  if (ncol(meth_data) != 7) stop("文件列数错误，预期7列，实际", ncol(meth_data), "列")
  
  # 转换关键列为数值型并检查异常值
  meth_data <- meth_data %>%
    mutate(
      # 强制转换base为整数（基因组位置必须为整数）
      base = as.integer(base),
      # 强制转换coverage为整数（覆盖度不能为非整数）
      coverage = as.integer(coverage),
      # 强制转换freqC/freqT为数值（甲基化频率）
      freqC = as.numeric(freqC),
      freqT = as.numeric(freqT)
    ) %>%
    # 过滤无效行（base、coverage、freqC有NA的行）
    filter(
      !is.na(base),        # 排除无效基因组位置
      !is.na(coverage),    # 排除无效覆盖度
      !is.na(freqC),       # 排除无效甲基化频率
      coverage >= 0        # 覆盖度不能为负数
    )
  
  # 验证数据有效性（清洗后无有效行则报错）
  if (nrow(meth_data) == 0) stop("文件无有效数据行（可能存在格式错误）: ", file_path)
  
  # 3. 修正strand值（保持原逻辑）
  meth_data <- meth_data %>%
    mutate(
      strand = case_when(
        strand %in% c("F", "+") ~ "+",
        strand %in% c("R", "-") ~ "-",
        TRUE ~ "*"  # 未知链标记为*
      )
    )
  
  # 4. 计算染色体长度（严格确保为数值）
  current_chr <- unique(meth_data$chr)
  if (length(current_chr) != 1) stop("文件包含多个染色体，应为单染色体文件: ", paste(current_chr, collapse = ","))
  cat("当前染色体:", current_chr, "\n")
  
  # 优先从.fai获取长度，否则用数据中最大base
  if (!is.null(fai_path) && file.exists(fai_path)) {
    fai <- read.table(fai_path, col.names = c("chr", "length", "offset", "linebases", "linewidth"), stringsAsFactors = FALSE)
    chrom_lengths <- fai$length[fai$chr == current_chr]
    if (length(chrom_lengths) == 0 || is.na(chrom_lengths)) {
      cat(".fai中未找到当前染色体，改用数据最大位置作为长度\n")
      chrom_lengths <- max(meth_data$base, na.rm = TRUE)
    }
  } else {
    chrom_lengths <- max(meth_data$base, na.rm = TRUE)
  }
  
  # 检查染色体长度有效性
  if (is.na(chrom_lengths) || chrom_lengths <= 0) stop("染色体长度无效（非数值或<=0）: ", chrom_lengths)
  cat("染色体长度:", chrom_lengths, "bp\n")
  
  # 5. 创建2kb窗口（确保start和end为数值）
  window_starts <- seq(1, chrom_lengths, by = 2000)
  window_ends <- pmin(seq(2000, chrom_lengths + 1999, by = 2000), chrom_lengths)
  # 检查窗口坐标是否为数值
  if (!is.numeric(window_starts) || !is.numeric(window_ends)) stop("窗口坐标非数值，无法创建窗口")
  
  windows <- GRanges(
    seqnames = current_chr,
    ranges = IRanges(start = window_starts, end = window_ends)
  )
  cat("生成窗口数:", length(windows), "\n")
  
  # 6. 计算窗口甲基化水平（保持原逻辑）
  gr <- GRanges(
    seqnames = meth_data$chr,
    ranges = IRanges(start = meth_data$base, width = 1),
    strand = meth_data$strand,
    coverage = meth_data$coverage,
    methylated = (meth_data$freqC / 100) * meth_data$coverage  # 直接计算甲基化位点数量
  )
  
  overlaps <- findOverlaps(windows, gr)
  window_df <- data.frame(
    chr = as.character(seqnames(windows)),
    start = start(windows),
    end = end(windows),
    context = context,
    total_coverage = 0,
    avg_methylation = NA_real_
  )
  
  for (i in seq_along(windows)) {
    idx <- subjectHits(overlaps[queryHits(overlaps) == i])
    if (length(idx) > 0) {
      total_meth <- sum(gr$methylated[idx])
      total_cov <- sum(gr$coverage[idx])
      window_df$total_coverage[i] <- total_cov
      window_df$avg_methylation[i] <- ifelse(total_cov > 0, (total_meth / total_cov) * 100, NA)
    }
  }
  return(window_df)
}

# ======================
# 主执行流程（修改为你的文件路径）
# ======================
setwd("/share/org/YZWL/yzwl_liuyy/project/lwm_test/rice_methy_analyze/02_bismark.2/04_bismark_methylation_hap2/windows/YZ001/1")

# 输入文件（确保是按染色体拆分后的文件）
files <- c(
  CG = "YZ001_1_clean_R1_bismark_bt2_pe.deduplicated.CX_report.txt_CG_methykit_Chr01.txt",
  CHG = "YZ001_1_clean_R1_bismark_bt2_pe.deduplicated.CX_report.txt_CHG_methykit_Chr01.txt",
  CHH = "YZ001_1_clean_R1_bismark_bt2_pe.deduplicated.CX_report.txt_CHH_methykit_Chr01.txt"
)

# .fai文件路径（若无可设为NULL）
fai_path <- "./YZ001_hap2.chr1.fai"  # 确保该文件中包含Chr01的长度

# 执行计算（添加错误捕获，方便定位问题文件）
results <- lapply(names(files), function(ctx) {
  cat("\n===== 处理类型:", ctx, "=====\n")
  tryCatch({
    calculate_methylation_windows(files[ctx], ctx, fai_path)
  }, error = function(e) {
    cat("处理", ctx, "时出错:", e$message, "\n")
    NULL  # 出错时返回NULL，不中断整体流程
  })
}) %>% bind_rows()

# 结果输出与可视化（保持原逻辑）
if (nrow(results) > 0) {
  cat("\n结果示例：\n")
  print(head(results, 3))
  print(tail(results, 3))
  
  write.csv(results, "rice_2kb.Chr1.methylation_full_coverage.csv", row.names = FALSE)
  cat("结果已保存至: rice_2kb.Chr1.methylation_full_coverage.csv\n")
  
  # 可视化（略，保持原逻辑）
} else {
  cat("无有效结果生成，请检查输入文件格式！\n")
}

