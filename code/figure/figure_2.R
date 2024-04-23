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
physiology_collection0 <- physiology_collection[!duplicated(physiology_collection$condition_unique),]
combine0 <- combine[, colnames(combine) %in% physiology_collection0$sampleID]

# input the compartment annotation
compartment_annotation <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/compare_compartment_annotation_with_and_without_manual_curation.xlsx")
compartment_annotation['fold'] <- compartment_annotation$annotation_curation / compartment_annotation$annotation_combine
# remove 'membrane'? as it actually covers "membrane protein" from different organelles
compartment_annotation <- compartment_annotation[compartment_annotation$compartment !="membrane",]
compartment0 <- compartment_annotation[compartment_annotation$compartment !="mitochondrion_unassigned",]
compartment0 <- compartment0[compartment0$compartment !="cytoplasm",]

compartment0 <- compartment0[order(compartment0[[4]], decreasing = TRUE), ]

# bar plot to analyze the protein num in each main component
ggplot(data=compartment0[1:10,], aes(x = reorder(compartment, annotation_curation), y = annotation_curation)) + 
  geom_bar(stat="identity",fill="#619CFF") +
  xlab("Cellular Component") + 
  ylab("Protein count") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_flip()

# bar plot to analyze the protein num in each compartment before and after curation based on experimental evidance
compartment0 <- compartment0[compartment0$fold != 1, ]


long_DF <- compartment0 %>% gather(type, count, 3:4)
long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, aes(fill=type, y=count, x=compartment)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Cellular Component") + 
  ylab("Protein count)") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_flip()


# calculate the organelle mass fraction variance 
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
sd0 <- apply(subset(ProMassRatio, select = 2:276), 1, sd, na.rm=TRUE) 
mean0 <- apply(subset(ProMassRatio, select = 2:276), 1, mean, na.rm=TRUE) 
# generate new dataframe
organelle_df = data.frame(compartment=ProMassRatio$compartment, mean=mean0, sd=sd0)
# check the main organelle
organelle <- c('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet', 'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'plasma membrane')
organelle_main <- organelle_df[organelle_df$compartment %in% organelle,]
# bar plot
ggplot(organelle_main) +
  geom_bar( aes(x= reorder(compartment, -mean), y=mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=compartment, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=0.5) +
  ylab("Mass fraction of protein") + 
  xlab("") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# calculate the correlation of organelle mass across unique conditions
physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")
physiology_collection0 <- physiology_collection[!duplicated(physiology_collection$condition_unique),]

ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]


ProMassRatio_ss <- ProMassRatio1[, colnames(ProMassRatio1) %in% physiology_collection0$sampleID]


# all correlation analysis of different samples
M <- cor(ProMassRatio_ss, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df1 <- data.frame(value=ss)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0.75, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Correlation coefficient for all component",
       y = "Density") 



# if based on the main organelle
ProMassRatio2 <- ProMassRatio1[ProMassRatio1$compartment %in% organelle, ]
ProMassRatio_ss <- ProMassRatio2[, colnames(ProMassRatio2) %in% physiology_collection0$sampleID]

# correlation analysis of different samples
M <- cor(ProMassRatio_ss, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df2 <- data.frame(value=ss)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0.75, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Correlation coefficient at main organelle level",
       y = "Density") 



# if based on the suborganelle
ProMassRatio3 <- ProMassRatio1[!(ProMassRatio$compartment %in% organelle), ]
ProMassRatio_ss <- ProMassRatio3[, colnames(ProMassRatio3) %in% physiology_collection0$sampleID]

# correlation analysis of different samples
M <- cor(ProMassRatio_ss, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df3 <- data.frame(value=ss)
# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0.75, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Correlation coefficient at sub-organelle level",
       y = "Density") 

# combine the above three result together
df1$type = "All component"
df2$type = "Main organelle"
df3$type = "Sub-organelle"

updated <- rbind(df1, df2, df3)
ggplot(updated, aes(x=value, color=category, fill=category)) +
  geom_density(alpha=0.3)

