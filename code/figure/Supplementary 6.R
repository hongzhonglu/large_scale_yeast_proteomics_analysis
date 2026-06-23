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



######################################################
# jianye dataset
# classify into two group
physiology_jianye <- physiology_collection[str_detect(physiology_collection$source,"sysbio_Jianye"),]
physiology_rosemary <- physiology_collection[str_detect(physiology_collection$condition_unique,"@NH4@N_limit@C_N_ratio=30"),]
physiology_rosemary <- physiology_rosemary[str_detect(physiology_rosemary$sampleID,"prot\\."),]

#physiology_Ibrahim <- physiology_collection[str_detect(physiology_collection$sampleID,"Chemostats_C_limit"),]
#physiology_yihui <- physiology_collection[str_detect(physiology_collection$sampleID,"sce_FY4_C"),]

# compartment
# ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/all_organelle_fraction_test.xlsx") # update on 2/10/2026

ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]



# group1 select
physiology <- physiology_jianye
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$S27_carbon_limit),]
Pro_mass_select0 <- t(Pro_mass_select[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select$compartment
Pro_mass_select0 <- as.data.frame(Pro_mass_select0)
Pro_mass_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)



# scatter plot
ggplot(Pro_mass_select0, aes(x=growth, y=ribosome)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Protein mass fraction in ribosome") +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(ribosome ~ growth, data = Pro_mass_select0)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6




# scatter plot
ggplot(Pro_mass_select0, aes(x=growth, y=nucleolus)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Protein mass fraction in nucleolus") +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(nucleolus ~ growth, data = Pro_mass_select0)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  








######################################################
# rosemary dataset
# group1 select
physiology <- physiology_rosemary
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$prot.21),]
Pro_mass_select0 <- t(Pro_mass_select[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select$compartment
Pro_mass_select0 <- as.data.frame(Pro_mass_select0)
Pro_mass_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)



# scatter plot
ggplot(Pro_mass_select0, aes(x=growth, y=ribosome)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Protein mass fraction in ribosome") +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(ribosome ~ growth, data = Pro_mass_select0)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  



# scatter plot
ggplot(Pro_mass_select0, aes(x=growth, y=nucleolus)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Protein mass fraction in nucleolus") +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(nucleolus ~ growth, data = Pro_mass_select0)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  




