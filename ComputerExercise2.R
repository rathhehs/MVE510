load("genome1.rdata")
load("genome2.rdata")
load("genome3.rdata")

genome1.subset=genome1[1:5000,]
ref.subset=reference[1:5000]

# 定义计算 p 值的函数
calculate_p_value <- function(y_i, N_i, p_error) {
  # 如果覆盖深度为零，返回 NA
  if (N_i == 0) {
    return(NA)
  }
  # 使用 binom.test 计算 p 值
  result <- binom.test(y_i, N_i, p_error, alternative = "greater")
  return(result$p.value)
}

genome2.length <- nrow(genome2)
mismatches <- vector(length = genome2.length)
coverage <- vector(length = genome2.length)
p_values <- vector(length = genome2.length)  # 用于存储 p 值

# 设置测序错误率 p_error（例如 Illumina 数据的平均错误率）
p_error <- 0.01

# 计算每个位点的错配数和覆盖深度，并计算 p 值
for (pos in 1:genome2.length) {
  ref_base <- reference[pos]
  
  genome_bases <- names(genome2[pos, -1])[genome2[pos, -1] > 0]
  genome_counts <- genome2[pos, -1][genome1[pos, -1] > 0]
  
  coverage[pos] <- sum(genome_counts)
  
  mismatches[pos] <- sum(genome_counts[genome_bases != ref_base])
  
  # 计算 p 值
  p_values[pos] <- calculate_p_value(mismatches[pos], coverage[pos], p_error)
}

# 标记显著的突变位置（例如使用 0.05 的显著性水平）
significant_positions <- which(p_values < 0.05)

# 创建一个包含显著位置和对应 p 值的数据框
significant_data <- data.frame(
  Position = significant_positions,
  P_value = p_values[significant_positions]
)

# 保存到 CSV 文件
write.csv(significant_data, file = "genome2.csv", row.names = FALSE)






