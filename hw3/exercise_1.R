####  1.2  #### 

# Data reading
setwd("~/MEGA/Master/BHT/Assignments/hw3/data_for_part_1/")

df <- read.table("normalized_data.txt", h = F)
names(df) <- c("D1", "D2", "D3", "D4", "D5", "C1", "C2", "C3", "C4","C5")


# As there is no evidence of that we can safely reject a normal distribution for all genes in both 
# groups, we will use a t.test to find whereas there are significant differences in gene expression
my_pvals <- apply(df, 1, function (x) t.test(x[1:5], x[6:10])$p.value)


####  1.3  #### 

# Finding the genes that are differentially expressing
diff_genes <- length(which(my_pvals < 0.05)) # 1911
diff_genes <- sum(my_pvals < 0.05)

# Expected FP = threshold% used in the t.test, in this case, 0.05% of the t. tests that reject the H0
# are wrong 
FP1 <- round(0.05 * diff_genes, 0) # 96


####  1.4  #### 

# Bonferroni multiple testing correction
my_bonferroni_pvals <- p.adjust(my_pvals, method = "bonferroni")

# Number of genes with p-value < 0.2 with Bonferroni correction
sum(my_bonferroni_pvals < 0.2) # 0

# BH multiple testing correction
my_BH_pvals <- p.adjust(my_pvals, method = "BH")

# Number of genes with p-value < 0.2 with BH correction
sum(my_BH_pvals < 0.2) # 12

# Expected FDR = threshold% used in the t.test, in this case, 0.2% of the t. tests
# are expected to be FP
FP2 <- round(0.2 * sum(my_BH_pvals < 0.2), 0) # 2


####  1.5  #### 

# Addition of means
df$D_mean <- rowMeans(df[, c(1:5)])
df$C_mean <- rowMeans(df[, c(6:10)])

# Calculate the log2 foldchange for each gene using this formula:
# foldchange = log2(mean(hiv)) - log2(mean(control)) # 
my_foldchanges <- log2(df[, 11]) - log2(df[, 12])


####  1.6  #### 

# Report the fold changes for the genes with a FDR < 0.2.
# Are there most up or down-regulated genes in the HIV subset?
which(my_BH_pvals < 0.2)

genes_FDR_lower_0.2 <- log2(df[which(my_BH_pvals < 0.2), 11]) - 
  log2(df[which(my_BH_pvals < 0.2), 12])

sum(genes_FDR_lower_0.2 > 0) # 12 genes upregulated

# It is interesting to point out that the 12 genes with a significative fold change are upregulated 
# in the HIV patients, where the gene with most fold change has a 0.6 fold change, this is, it is 50%
# more expressed than the same gene in the control condition. However, the gene with the lower 
# fold change has 0.04.

