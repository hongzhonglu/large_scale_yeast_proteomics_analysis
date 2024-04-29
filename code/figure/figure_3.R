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
physiology_rosemary <- physiology_collection[str_detect(physiology_collection$condition_unique,"@NH4@N_limit@C_N_ratio=30"),]
physiology_rosemary <- physiology_rosemary[str_detect(physiology_rosemary$sampleID,"prot\\."),]
physiology_Ibrahim <- physiology_collection[str_detect(physiology_collection$sampleID,"Chemostats_C_limit"),]

physiology_yihui <- physiology_collection[str_detect(physiology_collection$sampleID,"sce_FY4_C"),]

# compartment
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]



# group1 select
physiology <- physiology_jianye
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select0 <- t(Pro_mass_select[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select$compartment
Pro_mass_select0 <- as.data.frame(Pro_mass_select0)
Pro_mass_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)



# scatter plot
ggplot(Pro_mass_select0, mapping = aes(x=growth, y=mitochondrion)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth",
       y = "mitochondrion") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold"))

# plot for the main organelle
organelle <- c('mitochondrion', 'nucleus', 'cytosol', 'ribosome')
df_c <- Pro_mass_select0[, c("growth", organelle)]
long_DF <- df_c %>% gather(type, mass_fraction, 2:5)
long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))



#plot for the organelle with smaller fraction
organelle <- c('lipid droplet', 'fungal-type vacuole', 'endoplasmic reticulum', 'Golgi apparatus','peroxisome','endosome')
df_c <- Pro_mass_select0[, c("growth", organelle)]
long_DF <- df_c %>% gather(type, mass_fraction, 2:(1+length(organelle)))
long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))





# calculate the correlation of all component and growth rate
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$S27_carbon_limit),]
Pro_mass_select <- Pro_mass_select[Pro_mass_select$S27_carbon_limit >0.01,]
Pro_mass_select0 <- t(Pro_mass_select[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select$compartment
Pro_mass_select0 <- as.data.frame(Pro_mass_select0)
Pro_mass_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)
r <- cor(Pro_mass_select0)
r_growth <- r["growth",]
r_df <- as.data.frame(r_growth)
colnames(r_df) <- "cor"
r_df$component <- rownames(r_df)
r_df <- r_df[!is.na(r_df$cor), ]
r_df <- r_df[r_df$component !="growth", ]

#bar plot
ggplot(r_df, aes(x=reorder(component, cor), y=cor, fill=component)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Cellular Component") + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position="none")+
  coord_flip()
