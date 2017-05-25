####  1.2  #### 

# Data reading
setwd("C:/Users/Juanma/Desktop/Block 4/HIGH TROUGHPUT/HW/HW_3/data_for_hw3/data_for_part_1")
setwd("~/HT/HW/HW3/data_for_hw3/data_for_part_1/")

df <- read.table("normalized_data.txt", h = F)
names(df) <- c("D1", "D2", "D3", "D4", "D5", "C1", "C2", "C3", "C4","C5")


# Addition of means
df$D_mean <- rowMeans(df[, c(1:5)])
df$C_mean <- rowMeans(df[, c(6:10)])

head(df)

# Check for normality
par(mfrow = c(3, 4))
apply(df[, c(1:10)], 2, function (x) hist(x))

# Assuming normal distributions for all patient groups, perform a t.test
my_pvals <- apply(df[, c(11, 12)], 1, function (x) t.test(x)$p.value)

# Finding the genes that are NOT differentially expressing
which(my_pvals > 0.05) # 40


####  1.3  #### 

# Expected FDR = threshold% used in the t.test, in this case, 0.05% of the t. tests
# are expected to be FP
my_FDR_1 <- round(0.05 * length(my_pvals), 0) # 1114

# Genes with a p-value < 0.05
my_sig_pvals <- length(which(my_pvals <= 0.05)) # 22243


####  1.4  #### 

# Bonferroni multiple testing correction
my_bonferroni_pvals <- p.adjust(my_pvals, method = "bonferroni")

# Number of genes with p-value < 0.2 with Bonferroni correction
length(which(my_bonferroni_pvals < 0.2)) # 37

# Holm multiple testing correction
my_holm_pvals <- p.adjust(my_pvals, method = "BH")

# Number of genes with p-value < 0.2 with Holm correction
length(which(my_holm_pvals < 0.2)) # 37

# Expected FDR = threshold% used in the t.test, in this case, 0.2% of the t. tests
# are expected to be FP
my_FDR_2 <- round(0.2 * length(my_pvals), 0) # 4457


####  1.5  #### 

# Calculate the log2 foldchange for each gene using this formula:
            # foldchange = log2(mean(hiv)) - log2(mean(control)) # 
my_foldchanges <- log2(df[, 11]) - log2(df[, 12])


####  1.6  #### 

# Report the fold changes for the genes with a FDR < 0.2 (calculated with Bonferroni
# correction). Are there most up or down-regulated genes in the HIV subset?
which(my_bonferroni_pvals < 0.2)

genes_FDR_lower_0.2 <- log2(df[which(my_bonferroni_pvals < 0.2), 11]) - 
                       log2(df[which(my_bonferroni_pvals < 0.2), 12])

length(which(genes_FDR_lower_0.2 > 0)) # 18 genes upregulated
length(which(genes_FDR_lower_0.2 < 0)) # 19 genes downregulated


