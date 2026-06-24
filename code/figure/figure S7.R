library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)

# part 1 proteomics analysis-main part
# data and sample
combine <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/mass_fraction_combine.xlsx")


# analyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine[,2:276]))
# If you want the result as a dataframe
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine$gene)
gene_remove <- rownames(na_counts_df)[which(na_counts_df$row_NA_count >=275)]
combine <- combine[!(combine$gene %in% gene_remove),]


physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")


# extract the chemostat with different growth rate
physiology_collection <- physiology_collection[!is.na(physiology_collection$growth_mode),]
physiology_collection1 <- physiology_collection[str_detect(physiology_collection$growth_mode,"chemostat"), ]
physiology_collection1 <- physiology_collection1[!str_detect(physiology_collection1$Strain,"orientalis"), ]

# classify into two group
physiology_jianye <- physiology_collection[str_detect(physiology_collection$source,"sysbio_Jianye"),]



### randomly selected proteins and check the related tendencies
# compartment
physiology <- physiology_jianye
combine_select <-  combine[, colnames(combine) %in% c("gene",physiology$sampleID)]

combine_select0 <- t(combine_select[,-1])
colnames(combine_select0) <- combine_select$gene
combine_select0 <- as.data.frame(combine_select0)
combine_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)
# remove the column with all NA values
combine_select0 <- combine_select0[, colSums(is.na(combine_select0)) < nrow(combine_select0)]
gene_list <- colnames(combine_select0)
gene_list <- gene_list[which(gene_list !='growth')]


# randomly selected genes
ss <- 1000
cycle <- seq(from = 1, to = ss, by = 1)
new_df <- data.frame(growth=combine_select0$growth)

for (x in cycle){
  print(x)
  col0 <- paste("mass_fraction", x)
  number_gene <- ss
  gene_select  <- sample(gene_list, number_gene)
  combine_subset <- combine_select0[, colnames(combine_select0) %in%gene_select]
  new <- rowSums(combine_subset[,1:number_gene], na.rm = TRUE)
  new_df[,col0] <- new
}

new_df1 <- new_df[, !colnames(new_df) %in% c('growth')]
# calculate the average value
new11 <- rowSums(new_df1[,1:length(cycle)], na.rm = TRUE)/length(cycle)

SD <- apply(new_df1[,1:length(cycle)],1,sd)

final_df <- data.frame(growth=combine_select0$growth, mass_fraction=new11, sd=SD)


ggplot(final_df, mapping = aes(x=growth, y=mass_fraction)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction for 1000 proteins") +
  theme_bw() +
  geom_smooth() +
  ylim(0,0.4) +
  geom_errorbar(aes(ymin=mass_fraction - sd, ymax=mass_fraction + sd), width=.01, 
                position=position_dodge(0.05)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))





# randomly selected genes
ss <- 500
cycle <- seq(from = 1, to = ss, by = 1)
new_df <- data.frame(growth=combine_select0$growth)

for (x in cycle){
  print(x)
  col0 <- paste("mass_fraction", x)
  number_gene <- ss
  gene_select  <- sample(gene_list, number_gene)
  combine_subset <- combine_select0[, colnames(combine_select0) %in%gene_select]
  new <- rowSums(combine_subset[,1:number_gene], na.rm = TRUE)
  new_df[,col0] <- new
}

new_df1 <- new_df[, !colnames(new_df) %in% c('growth')]
# calculate the average value
new11 <- rowSums(new_df1[,1:length(cycle)], na.rm = TRUE)/length(cycle)

SD <- apply(new_df1[,1:length(cycle)],1,sd)

final_df <- data.frame(growth=combine_select0$growth, mass_fraction=new11, sd=SD)


ggplot(final_df, mapping = aes(x=growth, y=mass_fraction)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction for 500 proteins") +
  theme_bw() +
  geom_smooth() +
  ylim(0,0.4) +
  geom_errorbar(aes(ymin=mass_fraction - sd, ymax=mass_fraction + sd), width=.01, 
                position=position_dodge(0.05)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))





# randomly selected genes
ss <- 200
cycle <- seq(from = 1, to = ss, by = 1)
new_df <- data.frame(growth=combine_select0$growth)

for (x in cycle){
  print(x)
  col0 <- paste("mass_fraction", x)
  number_gene <- ss
  gene_select  <- sample(gene_list, number_gene)
  combine_subset <- combine_select0[, colnames(combine_select0) %in%gene_select]
  new <- rowSums(combine_subset[,1:number_gene], na.rm = TRUE)
  new_df[,col0] <- new
}

new_df1 <- new_df[, !colnames(new_df) %in% c('growth')]
# calculate the average value
new11 <- rowSums(new_df1[,1:length(cycle)], na.rm = TRUE)/length(cycle)

SD <- apply(new_df1[,1:length(cycle)],1,sd)

final_df <- data.frame(growth=combine_select0$growth, mass_fraction=new11, sd=SD)


ggplot(final_df, mapping = aes(x=growth, y=mass_fraction)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction for 200 proteins") +
  theme_bw() +
  geom_smooth() +
  ylim(0,0.4) +
  geom_errorbar(aes(ymin=mass_fraction - sd, ymax=mass_fraction + sd), width=.01, 
                position=position_dodge(0.05)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
