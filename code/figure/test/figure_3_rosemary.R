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

# compartment
# ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/all_organelle_fraction_test.xlsx") # update on 2/10/2026


ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]



# group1 select
physiology <- physiology_rosemary
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$prot.21),]
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



# plot for the main organelle membrane
main_membrane <- c("mitochondrial outer membrane","mitochondrial inner membrane", "endoplasmic reticulum membrane",
                   "plasma membrane" )
df_m <- Pro_mass_select0[, colnames(Pro_mass_select0) %in% main_membrane ]
df_m$growth <- Pro_mass_select0$growth
long_DF <- df_m %>% gather(type, mass_fraction, 1:4)
long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))



# plot for the main organelle - mt
df_mt <- Pro_mass_select0[, str_detect(colnames(Pro_mass_select0), "^mitochondrial ")]
df_mt$growth <- Pro_mass_select0$growth
long_DF <- df_mt %>% gather(type, mass_fraction, 1:8)
long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),legend.title=element_text(size=12), legend.text=element_text(size=12))



# plot for the main organelle - nuclear
df_o <- Pro_mass_select0[, str_detect(colnames(Pro_mass_select0), "^nucl")]
col_select <- c("nucleosome","nucleolus", "nucleoplasm", "nuclear periphery", "nuclear chromosome", "nuclear envelope")
df_o <- df_o[, colnames(df_o) %in% col_select]
df_o$growth <- Pro_mass_select0$growth
#long_DF <- df_o %>% gather(type, mass_fraction, 1:6)
long_DF <- df_o %>%
  pivot_longer(
    cols = -growth,                    # all columns except growth
    names_to = "type",
    values_to = "mass_fraction"
  )

long_DF$type <- as.factor(long_DF$type)
ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),legend.title=element_text(size=12), legend.text=element_text(size=12))





















# plot for the main organelle membrane
# change all mass fraction into the relative value
#main_membrane <- c("fungal-type vacuole membrane" , "mitochondrial outer membrane","mitochondrial inner membrane", "endoplasmic reticulum membrane",
#                   "plasma membrane","nuclear inner membrane", "peroxisomal membrane" ,"nuclear outer membrane","Golgi membrane",
#                   "prospore membrane", "cellular bud membrane", "endosome membrane" )
main_membrane <- c("fungal-type vacuole membrane" , "mitochondrial outer membrane","mitochondrial inner membrane", "endoplasmic reticulum membrane",
                   "plasma membrane","nuclear inner membrane", "peroxisomal membrane" ,"nuclear outer membrane","Golgi membrane",
                   "endosome membrane" )
df_m <- Pro_mass_select0[, colnames(Pro_mass_select0) %in% main_membrane ]
df_m_ref = df_m[c(1:3),]
avearge_value <- apply(df_m_ref, 2, mean)

df_m_relative <- mapply('/', df_m, avearge_value)

df_m_relative <- as.data.frame(df_m_relative)
df_m_relative$growth <- Pro_mass_select0$growth
long_DF <- df_m_relative %>% gather(type, mass_fraction, 1:10)
long_DF$type <- as.factor(long_DF$type)

ggplot(long_DF, mapping = aes(x=growth, y=mass_fraction, colour=type)) +
  geom_point() + # geom_point(alpha = 2/10) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Growth rate (/h)",
       y = "Scaled mass fraction") +
  theme_bw() +
  geom_smooth()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),legend.title=element_text(size=12), legend.text=element_text(size=12)) +
  geom_hline(yintercept = 1, col = "black")





# calculate the correlation of all component and growth rate
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$prot.21),]
Pro_mass_select_F <- Pro_mass_select[Pro_mass_select$prot.21 >0.01,]
Pro_mass_select0 <- t(Pro_mass_select_F[,-1])
colnames(Pro_mass_select0) <- Pro_mass_select_F$compartment
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








# focus on one specific organelle-for example nucleolus
compartment_sce_curation <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/compartment_sce_curation.xlsx")

compartment_one <- compartment_sce_curation[compartment_sce_curation$compartment=="nucleolus",]
compartment_one <- compartment_sce_curation[compartment_sce_curation$compartment=="peroxisome",]


combine_one <- combine[combine$gene %in% compartment_one$gene,]
combine_one_select <-  combine_one[, colnames(combine_one) %in% c("gene",physiology$sampleID)]

# remove too much na in each row
na_counts_per_row <- rowSums(is.na(combine_one_select[,2:19]))
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine_one_select$gene)
gene_remove <- rownames(na_counts_df)[which(na_counts_df$row_NA_count >=18)]
combine_one_select <- combine_one_select[!(combine_one_select$gene %in% gene_remove),]

combine_one_select0 <- t(combine_one_select[,-1])
colnames(combine_one_select0) <- combine_one_select$gene
combine_one_select0 <- as.data.frame(combine_one_select0)
combine_one_select0$growth <- as.numeric(physiology$`dilution rate (/h)`)

r <- cor(combine_one_select0, use='pairwise.complete.obs') # also calculate the columns with NA
r_growth <- r["growth",]
r_df <- as.data.frame(r_growth)
colnames(r_df) <- "cor"
r_df$gene <- rownames(r_df)
r_df <- r_df[!is.na(r_df$cor), ]
r_df <- r_df[r_df$gene !="growth", ]
#bar plot
ggplot(r_df, aes(x=reorder(gene, -cor), y=cor, fill=gene)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Gene") + 
  theme(axis.text = element_text(size = 3), axis.title = element_text(size = 12))+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position="none")

