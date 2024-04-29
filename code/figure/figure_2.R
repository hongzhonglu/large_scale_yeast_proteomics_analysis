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



# if based on the main organelle
ProMassRatio2 <- ProMassRatio1[ProMassRatio1$compartment %in% organelle, ]
ProMassRatio_ss <- ProMassRatio2[, colnames(ProMassRatio2) %in% physiology_collection0$sampleID]
# correlation analysis of different samples
M <- cor(ProMassRatio_ss, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df2 <- data.frame(value=ss)


# if based on the suborganelle
ProMassRatio3 <- ProMassRatio1[!(ProMassRatio1$compartment %in% organelle), ]
ProMassRatio_ss <- ProMassRatio3[, colnames(ProMassRatio3) %in% physiology_collection0$sampleID]
# correlation analysis of different samples
M <- cor(ProMassRatio_ss, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df3 <- data.frame(value=ss)


# combine the above three result together
df1$type = "All component"
df2$type = "Main organelle"
df3$type = "Sub-organelle"

updated <- rbind(df1, df2, df3)
ggplot(updated, aes(x=value, color=type, fill=type)) +
  geom_density(alpha=0.3) +
  xlim(0.5, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "Correlation coefficient between samples",
       y = "Density") 

# heatmap of mass fraction of main organelle in each unique condition?
library("pheatmap")
#organelle_s <- c('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum', 'fungal-type vacuole','peroxisome','ribosome')
ProMassRatio2 <- ProMassRatio1[ProMassRatio1$compartment %in% organelle, ]
rownames0 <- ProMassRatio2$compartment
DF <- ProMassRatio2[,-1]
rownames(DF) <- rownames0
pheatmap(DF, scale="none",
         cutree_rows = 4,
         show_colnames =FALSE)

pheatmap(DF, scale="column",
         cutree_rows = 4,
         show_colnames =FALSE)


suorganelle <- c('mitochondrial outer membrane','mitochondrial inner membrane','mitochondrial intermembrane space','mitochondrial matrix')
ProMassRatio_suborganelle <- ProMassRatio1[ProMassRatio1$compartment %in% suorganelle, ]

rownames0 <- ProMassRatio_suborganelle$compartment
DF <- ProMassRatio_suborganelle[,-1]
rownames(DF) <- rownames0
pheatmap(DF, scale="none",
         cutree_rows = 4,
         show_colnames =FALSE)













# supplementary figure
# calculate the protein correlation under different samples
combine00 <- combine[,  colnames(combine) %in% c("gene",physiology_collection0$sampleID)]
combine00 <- combine[, str_detect(colnames(combine), "prot\\.")]
combine00$gene <- combine$gene
# analyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine00[, 1:42]))
# If you want the result as a dataframe
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine00$gene)
gene_remove <- rownames(na_counts_df)[which(na_counts_df$row_NA_count >=42)]
combine01 <- combine00[!(combine00$gene %in%gene_remove), ]
combine02 <- t(combine01[, c(1:42)])
colnames(combine02) <- combine01$gene
# all correlation analysis of different samples
M <- cor(combine02, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
ss <- as.vector(M[upper.tri(M)])
df_all <- data.frame(cor0=ss, type="all")
# input the compartment annotation
compartment_sce <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/compartment_sce_curation.xlsx")
compartment_sce <- compartment_sce[compartment_sce$compartment !="membrane",]
compartment_sce <- compartment_sce[compartment_sce$compartment !="mitochondrion_unassigned",]
compartment_sce <- compartment_sce[compartment_sce$compartment !="cytoplasm",]
compartment_sce <- compartment_sce[, c(2:3)]


organelle <- c('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet', 'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'plasma membrane','mitochondrial outer membrane','mitochondrial inner membrane', 'nucleolus','mitochondrial intermembrane space','mitochondrial matrix')
for (x in unique(organelle)){
  print(x)
  cc=x
  gene_m <- compartment_sce[compartment_sce$compartment ==cc,]
  gene_calcualted <- colnames(M)
  gene_m_selected <- gene_m$gene[which(gene_m$gene %in% gene_calcualted)]
  M_m <- M[gene_m_selected, gene_m_selected]
  ss_m <- as.vector(M_m[upper.tri(M_m)])
  df_m <- data.frame(cor0=ss_m, type=cc)
  
  
  df_random_sample  <- data.frame(cor0=sample(ss, length(ss_m)), type="random")
  df_combine_test <- rbind(df_all, df_m, df_random_sample)
  
  ggplot(df_combine_test, aes(x=cor0, color=type, fill=type)) +
    geom_density(alpha=0.3) +
    xlim(-1, 1) +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    geom_density(alpha = 0.5)+
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
    labs(x = "Correlation coefficient of protein pair",
         y = "Density") 
  ggsave(out <- paste('result/',x,'.png', sep = ""), width=8, height=6, dpi=600)
}

















# supplementary figure
# correlation between mass and protein volume
sce_protein_MW_and_volume <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/sce_protein_MW_and_volume.xlsx")
# scatter plot
ggplot(sce_protein_MW_and_volume, mapping = aes(x=MW,y=Total_Volume)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Protein molecular weight",
       y = "Protein 3D structure volume") +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold"))


# compare Volume ration and mass ratio
ProVolumeRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/volume_size_ratio_across_compartment_combine.xlsx")
ProVolumeRatio <- ProVolumeRatio[,2:277]
ProVolumeRatio[ProVolumeRatio <0.0000000000001] <- NA
ProVolumeRatio1 <- ProVolumeRatio[ProVolumeRatio$compartment !="cytoplasm", ]
ProVolumeRatio1 <- ProVolumeRatio1[ProVolumeRatio1$compartment !="mitochondrion_unassigned", ]

mass_vs_volume <- ProMassRatio1[,c(1,2)]
colnames(mass_vs_volume) <- c("compartment","mass_fraction")
mass_vs_volume$volume_fraction <- ProVolumeRatio1$`Glucose_phase_rep1(g/gDW)`

# scatter plot
ggplot(mass_vs_volume, aes(x=mass_fraction, y=volume_fraction)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Volume fraction of components") +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold"))

