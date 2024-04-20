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



# calculate the mass ratio of top 1000 proteins
ratio_all <- c()
for(i in 2:ncol(combine)){
  print(i)
  ss1 <- combine[,i]
  ss1 <- ss1[order(ss1[[1]], decreasing = TRUE), ]
  top1000 <- ss1[1:1000,]
  ratio1 <- sum(top1000[,1],na.rm = TRUE)/sum(ss1[,1], na.rm = TRUE)
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
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Mass ratio of top 1000 proteins",
       y = "Density") 


# analyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine[,2:276]))


# If you want the result as a dataframe
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine$gene)






# part2 protein compartment







# part3 protein properties - mass, volume, sectional area






