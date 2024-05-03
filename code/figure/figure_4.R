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
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$S27_carbon_limit),]
Pro_mass_select0 <- t(Pro_mass_select[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select$compartment
Pro_mass_select0 <- as.data.frame(Pro_mass_select0)
Pro_mass_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)
Pro_mass_select0$growth <- as.factor(Pro_mass_select0$growth)
# calculate the average value of each group
u.id <- unique(Pro_mass_select0$growth)
newdata <- t(sapply(unique(u.id),function(c.id){
  colMeans(Pro_mass_select0[Pro_mass_select0$growth == c.id,-157],na.rm = TRUE)
}))

rownames(newdata) <- u.id
newdata2 <- t(newdata)
M <- cor(newdata2, method = "pearson", use = "pairwise.complete.obs") # for each pair, only non-NA value was calculated
corrplot(M[1:9, 1:9], tl.pos='n', type = "upper", shade.lwd=0.0001)
corrplot.mixed(M[1:9, 1:9], tl.pos='n', shade.lwd=0.0001)

write.table(newdata, "newdata.txt", sep = "\t")






'''
data <- data.frame(name=as.factor(c("a","b","c","a","c")),
                   v1=c(1,3,2,6,3),
                   v2=c(2,8,5,0,9),
                   v3=c(7,6,0,6,4),
                   v4=c(9,4,1,2,7),
                   v5=c(3,8,9,1,5))

u.id <- unique(data$name)
newdata <- t(sapply(unique(u.id),function(c.id){
  colMeans(data[data$name == c.id,-1])
}))

rownames(newdata) <- u.id
'''