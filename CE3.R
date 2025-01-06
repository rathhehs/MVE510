# 3.1
x = read.table('counts_matrix.txt')
metadata = read.table('metadata.txt', sep='\t', header=TRUE)

countNonZero <- function(x) {
  counts <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    counts[i] <- sum(x[i, ] > 0)
  }
  counts
}

x.filtered <- x[counts >= 20, ]
proportion_filtered <- nrow(x.filtered) / nrow(x) > 0.25
nrow(x.filtered)

# 3.2
calculateCPM <- function(x) {
  x <- x + 1
  numCols <- ncol(x)
  numRows <- nrow(x)
  cpmMatrix <- matrix(0, numRows, numCols)
  
  for (col in seq_len(numCols)) {
    colSum <- sum(x[, col])
    for (row in seq_len(numRows)) {
      cpmMatrix[row, col] <- (x[row, col] / colSum) * 1e6
    }
  }
  
  return(cpmMatrix)
}

filteredCPM <- calculateCPM(x.filtered)
filteredLogCPM <- log(filteredCPM)

# 3.3
x.filtered.log <- log(1 + x.filtered)

boxplot(filteredLogCPM, col = "green", ylab = "Log-transformed CPM", add = FALSE)
title("Boxplot of all genes in each sample with normalization")

boxplot(x.filtered.log, col = "green", ylab = "Log-transformed x.filtered", add = FALSE)
title("Boxplot of all genes in each sample without normalization")

plot(filteredLogCPM[, 1], filteredLogCPM[, 41], col = "green", 
     main = "Scatter plot for Sample 1 and Sample 41", 
     xlab = "Sample 1", ylab = "Sample 41")

# 3.4
gene1_cpm_values <- t(filteredCPM[1,])
gene1_logcpm_values <- t(filteredLogCPM[1,])

boxplot(gene1_cpm_values[1:40], gene1_cpm_values[41:80], names = c("Not IBD", "CD"))
title("Boxplot of Gene 1 with CPM Data")

boxplot(gene1_logcpm_values[1:40], gene1_logcpm_values[41:80], names = c("Not IBD", "CD"))
title("Boxplot of Gene 1 with LogCPM Data")

diagnosis_factor <- as.factor(metadata[,"diagnosis"])
diagnosis_factor <- relevel(diagnosis_factor, ref = "Not IBD")
sex_factor <- as.factor(metadata[,"Sex"])
age_at_diagnosis <- metadata[,"age.at.diagnosis"]

fit_diagnosis_model <- lm(gene1_logcpm_values[1,] ~ diagnosis_factor)
fit_full_model <- lm(gene1_logcpm_values[1,] ~ age_at_diagnosis + sex_factor + diagnosis_factor)

summary(fit_diagnosis_model)$coefficients
summary(fit_full_model)$coefficients

# 3.5
gene.length = nrow(filteredLogCPM)
fit1.pval = fit1.coeff = fit2.pval = fit2.coeff = vector(length = gene.length, mode = "double")

for (i in 1:gene.length) {
  fit1.i = lm(filteredLogCPM[i,] ~ diagnosis_factor)
  fit2.i = lm(filteredLogCPM[i,] ~ age_at_diagnosis + sex_factor + diagnosis_factor)
  
  fit1.pval[i] = summary(fit1.i)$coefficients["diagnosis_factorCD", "Pr(>|t|)"]
  fit1.coeff[i] = summary(fit1.i)$coefficients["diagnosis_factorCD", "Estimate"]
  
  fit2.pval[i] = summary(fit2.i)$coefficients["diagnosis_factorCD", "Pr(>|t|)"]
  fit2.coeff[i] = summary(fit2.i)$coefficients["diagnosis_factorCD", "Estimate"]
}

fit1.padjust = data.frame(p.adjust(fit1.pval, method = "fdr"), fit1.coeff)
fit2.padjust = data.frame(p.adjust(fit2.pval, method = "fdr"), fit2.coeff)

rownames(fit1.padjust) = rownames(fit2.padjust) = rownames(filteredLogCPM)
colnames(fit1.padjust) = colnames(fit2.padjust) = c("p-adjust", "coefficients")

cutoff = 0.05

fit1.diff.coeff = fit1.padjust[fit1.padjust[, 1] < cutoff, 2]
fit2.diff.coeff = fit2.padjust[fit2.padjust[, 1] < cutoff, 2]

num.fit1.diff = length(fit1.diff.coeff)
num.fit2.diff = length(fit2.diff.coeff)

fit1.coeff.up.num = sum(fit1.diff.coeff > 0)
fit1.coeff.down.num = num.fit1.diff - fit1.coeff.up.num
fit2.coeff.up.num = sum(fit2.diff.coeff > 0)
fit2.coeff.down.num = num.fit2.diff - sum(fit2.diff.coeff > 0)

fit1.most.gene = rownames(fit1.padjust)[which.min(fit1.padjust[, 1])]
fit2.most.gene = rownames(fit2.padjust)[which.min(fit2.padjust[, 1])]

coeff.most.gene.fit1 = fit1.padjust["ENSG00000185499", 2]
coeff.most.gene.fit2 = fit2.padjust["ENSG00000185499", 2]

log.fold.change = log(sum(filteredLogCPM["ENSG00000185499", 41:80]) / sum(filteredLogCPM["ENSG00000185499", 1:40]))
log.fold.change.full = log(sum(filteredCPM["ENSG00000185499", 41:80]) / sum(filteredCPM["ENSG00000185499", 1:40]))

fit1.padjust.sorted = fit1.padjust[order(fit1.padjust[, 1]), ]
fit2.padjust.sorted = fit2.padjust[order(fit2.padjust[, 1]), ]

rownames(fit1.padjust.sorted)[1:6]
fit1.padjust.sorted[1:6, 2] > 0

rownames(fit2.padjust.sorted)[1:6]
fit2.padjust.sorted[1:6, 2] > 0

fit2.age.pval = fit2.sex.pval = vector(length = gene.length, mode = "double")

for (i in 1:gene.length) {
  fit2.i = lm(filteredLogCPM[i, ] ~ age_at_diagnosis + sex_factor + diagnosis_factor)
  fit2.age.pval[i] = summary(fit2.i)$coefficients["age_at_diagnosis", "Pr(>|t|)"]
  fit2.sex.pval[i] = summary(fit2.i)$coefficients["sex_factorMale", "Pr(>|t|)"]
}

fit2.age.padjust = p.adjust(fit2.age.pval)
fit2.sex.padjust = p.adjust(fit2.sex.pval)

rownames(filteredLogCPM)[fit2.sex.padjust == min(fit2.sex.padjust)]

# 3.6
library(gplots)
library(RColorBrewer)

xsig = as.matrix(filteredLogCPM[rownames(fit1.padjust.sorted)[1:100],])
mycols = rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))
column.cols = c("purple", "orange")[metadata$diagnosis]

heatmap.2(xsig, trace='none', col=mycols, main='The 100 most significant genes', ColSideColors=column.cols)

# 3.7
pca = prcomp(t(filteredLogCPM))
summary(pca)

PC1 = pca$x[, 1]
PC2 = pca$x[, 2]
PC3 = pca$x[, 3]

colors1 = c("red", "blue")[factor(metadata$diagnosis, levels = c("CD", "Not IBD"))]

plot(PC1, PC2, 
     main = "Scatter plot of PC1 and PC2", 
     xlab = "PC1", ylab = "PC2", 
     col = colors1, pch = 19, cex = 1.5)

plot(PC1, PC3, 
     main = "Scatter plot of PC1 and PC3", 
     xlab = "PC1", ylab = "PC3", 
     col = colors1, pch = 19, cex = 1.5)

plot(PC2, PC3, 
     main = "Scatter plot of PC2 and PC3", 
     xlab = "PC2", ylab = "PC3", 
     col = colors1, pch = 19, cex = 1.5)

gender.cols = c("red", "blue")[factor(metadata$Sex, levels = c("Male", "Female"))]

plot(PC2, PC3, 
     main = "Scatter plot of PC2 and PC3", 
     xlab = "PC2", ylab = "PC3", 
     col = gender.cols, pch = 19, cex = 1.5)

plot(PC1, PC3, 
     main = "Scatter plot of PC1 and PC3", 
     xlab = "PC1", ylab = "PC3", 
     col = gender.cols, pch = 19, cex = 1.5)

plot(PC1, PC2, 
     main = "Scatter plot of PC1 and PC2", 
     xlab = "PC1", ylab = "PC2", 
     col = gender.cols, pch = 19, cex = 1.5)





