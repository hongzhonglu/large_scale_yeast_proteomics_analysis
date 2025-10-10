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


# supplementary figure 5
# correlation between mass and protein volume
sce_protein_MW_and_volume <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/sce_protein_MW_and_volume.xlsx")
# scatter plot
ggplot(sce_protein_MW_and_volume, mapping = aes(x=MW,y=Total_Volume)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Protein molecular weight (MW)",
       y = "Protein 3D structure volume (nm^3)") +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))


# 拟合线性模型
model <- lm(Total_Volume ~ MW, data = sce_protein_MW_and_volume)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6






# get the protein abundance under batch condition
library(dplyr)
batch_glucose <- combine %>% select(1:4)
batch_glucose$row_mean <- rowMeans(batch_glucose[, 2:4])

sce_protein_MW_and_volume$mass_fraction <- getSingleReactionFormula(batch_glucose$row_mean, batch_glucose$gene, sce_protein_MW_and_volume$locus)


# 删除 Age 列含 NA 的行
df_clean <- sce_protein_MW_and_volume[!(sce_protein_MW_and_volume$mass_fraction=="NA"), ]
df_clean$mass_fraction <- as.numeric(df_clean$mass_fraction)
df_clean$protein_relative_copy <- df_clean$mass_fraction / df_clean$MW
df_clean$protein_relative_volume <- df_clean$protein_relative_copy * df_clean$Total_Volume
sum(df_clean$mass_fraction)
sum(df_clean$protein_relative_volume)
df_clean$volume_fraction <- df_clean$protein_relative_volume/sum(df_clean$protein_relative_volume)

# scatter plot
ggplot(df_clean, mapping = aes(x=mass_fraction,y=volume_fraction)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction (%)",
       y = "Volume fraction (%)") +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))


# 拟合线性模型
model <- lm(volume_fraction ~ mass_fraction, data = df_clean)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6









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


