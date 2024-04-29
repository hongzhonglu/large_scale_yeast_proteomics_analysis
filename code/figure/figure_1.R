library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)
# corrplot
# part 1 proteomics analysis-main part
# data and sample
library(readxl)
combine <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/mass_fraction_combine.xlsx")


# analyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine[,2:276]))
# If you want the result as a dataframe
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine$gene)
gene_remove <- rownames(na_counts_df)[which(na_counts_df$row_NA_count >=275)]

combine <- combine[!(combine$gene %in% gene_remove),]



physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")
physiology_collection0 <- physiology_collection[!duplicated(physiology_collection$condition_unique),]

combine0 <- combine[, colnames(combine) %in% physiology_collection0$sampleID]


# correlation analysis of different samples
M <- cor(combine0, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
corrplot(M[1:50, 1:50], tl.pos='n', type = "upper", shade.lwd=0.0001)
corrplot.mixed(M[1:10, 1:10], tl.pos='n', shade.lwd=0.0001)
ss <- as.vector(M[upper.tri(M)])
df <- data.frame(value=ss)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Correlation coefficient",
       y = "Density") 



# cluter analysis
library(plotly)
library(Rtsne) # tSNE in an acronym for t-Distributed Neighbor Embedding is a statistical method that is mainly used to visualize high-dimensional data
library(umap) # umap is similar to tSNE, but more efficient

combine11 <- t(combine0)
combine11[is.na(combine11)] <- 0 # here NA value was replaced as 0

label11 <- as.factor(physiology_collection0$source[1:150])
iris.umap = umap(combine11, n_components = 2, random_state = 15) 

layout <- iris.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, label11) 

fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~label11, type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Source')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig 



# calculate the number of existing gene in each sample?
gene_num_per_sample <- c()
colnames00 <- colnames(combine)
for(i in 2:ncol(combine)){
  print(i)
  ss1 <- combine[,c(1, i)]
  ss1 <- ss1[order(ss1[[2]], decreasing = TRUE), ]
  ss2 <- ss1[!is.na(ss1[[2]]),]
  num00 <- length(ss2$gene)
  gene_num_per_sample <- c(gene_num_per_sample,  num00)
}

gene_num_df <- data.frame(sample=colnames00[2:length(colnames00)], gene_num=gene_num_per_sample)
# plot density plot
ggplot(gene_num_df, aes(gene_num)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 15, face = "bold")) +
  labs(x = "Gene number in each dataset",
       y = "Density") 






# calculate the mass ratio of top 1000 proteins
ratio_all <- c()
for(i in 2:ncol(combine)){
  print(i)
  ss1 <- combine[,c(1, i)]
  ss1 <- ss1[order(ss1[[2]], decreasing = TRUE), ]
  top1000 <- ss1[1:1000,]
  ratio1 <- sum(top1000[,2],na.rm = TRUE)/sum(ss1[,2], na.rm = TRUE)
  ratio_all <- c(ratio_all, ratio1)
}

df <- data.frame(value=ratio_all)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0.65, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 15, face = "bold")) +
  labs(x = "Mass ratio of top 1000 most abundant proteins",
       y = "Density") 



# reanalyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine[,2:276]))
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine$gene)
gene_existence <- na_counts_df
gene_existence$exist_count <- 275-gene_existence$row_NA_count

# 计算累计值
# 绘制累积分布图
ggplot(data.frame(x = gene_existence$exist_count), aes(x)) +
  stat_ecdf(geom = "step") +
  xlab("Existence num in all samples") +
  ylab("ECDF") +
  #theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 218, linetype = "dotted", # top 1500 existing proteins in all samples
    color = "red", linewidth = 1
  )

# plot density plot
ggplot(gene_existence, aes(exist_count)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 15, face = "bold")) +
  labs(x = "Gene occurance number across samples",
       y = "Density") 


# select the most existing protein
gene_rank <- gene_existence[order(gene_existence[[2]], decreasing = TRUE), ]
gene_rank0 <- gene_rank[c(1:1500),]

ratio_all2 <- c()
for(i in 2:ncol(combine)){
  print(i)
  ss1 <- combine[,c(1, i)]
  top1000_exist <- ss1[ss1$gene %in% rownames(gene_rank0), ]
  ratio1 <- sum(top1000_exist[,2],na.rm = TRUE)/sum(ss1[,2], na.rm = TRUE)
  ratio_all2 <- c(ratio_all2, ratio1)
}

df <- data.frame(value=ratio_all2)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0.4, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 15, face = "bold")) +
  labs(x = "Mass ratio of 1500 most frequent proteins",
       y = "Density") 


# analyze the gene with fewer occurrence num
gene_existence_l <- gene_existence[gene_existence$exist_count <= 10,]
gene_list <- rownames(gene_existence_l)
paste0(gene_list,collapse = ",")

# heatmap for some specific genes
library("pheatmap")
combine_1500 <- combine[combine$gene %in% rownames(gene_rank0[1:10, ]),]
combine_1500_s <- combine_1500[, colnames(combine_1500) %in% c('gene',physiology_collection0$sampleID)]
rownames0 <- combine_1500_s$gene
DF <- combine_1500_s[,-1]
rownames(DF) <- rownames0
pheatmap(DF, scale="column",
         show_colnames =FALSE,
         show_rownames =TRUE)


combine_1500 <- combine[combine$gene %in% rownames(gene_rank0[gene_rank0$exist_count>=274, ]),]
combine_1500_s <- combine_1500[, colnames(combine_1500) %in% c('gene',physiology_collection0$sampleID)]
rownames0 <- combine_1500_s$gene
DF <- combine_1500_s[,-1]
rownames(DF) <- rownames0
pheatmap(DF, scale="column",
         show_colnames =FALSE,
         show_rownames =TRUE)