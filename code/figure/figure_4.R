library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)


physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")


# classify into two group
physiology_rosemary2 <- physiology_collection[str_detect(physiology_collection$source, "sysbio_Rosemary2"),]
physiology_yirong <- physiology_collection[str_detect(physiology_collection$condition_unique, "Copy number WT"),]

# compartment
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]



# group1 select
physiology <- physiology_yirong
Pro_mass_select <-  ProMassRatio1[, colnames(ProMassRatio1) %in% c("compartment",physiology$sampleID)]

#combine the same condition
condition <- str_trim(colnames(Pro_mass_select), side = "both")
condition <- str_replace_all(condition, "hr [:digit:]$", "")
colnames(Pro_mass_select) <- condition
Pro_mass_select <- Pro_mass_select[!is.na(Pro_mass_select$`Copy number WT 0 `),]
Pro_mass_select0 <- Pro_mass_select[, !colnames(Pro_mass_select) %in% c("compartment")]
# calculate the average of columns
df <- as.data.frame(sapply(split.default(Pro_mass_select0, names(Pro_mass_select0)), rowMeans))


df$compartment <- Pro_mass_select$compartment
colnames(df) <- c("WT_0",  "WT_12", "WT_16", "WT_4",  "WT_8",  "compartment"  ) 


ggplot(df, aes(x=WT_0 , y=WT_4)) +
  geom_point(size=2, shape=23) +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 4h") 

ggplot(df, aes(x=WT_0 , y=WT_8)) +
  geom_point(size=2, shape=23) +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 8h") 


ggplot(df, aes(x=WT_0 , y=WT_12)) +
  geom_point(size=2, shape=23) +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 12h") 



df$WT_8_vs_WT_0 <- (df$WT_8-df$WT_0)/df$WT_0
df$WT_8_vs_WT_0_abs <- abs(df$WT_8_vs_WT_0)

df_filter_mass <- df[df$WT_0 >=0.01, ]
df_filter_mass <- df_filter_mass[df_filter_mass$WT_8_vs_WT_0_abs >=0.2, ] 

ggplot(df_filter_mass, aes(x=reorder(compartment, WT_8_vs_WT_0), y=WT_8_vs_WT_0)) + 
  geom_bar(position="dodge", stat="identity", fill = "blue") +
  xlab("Cellular Component") + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position="none")+
  coord_flip()
  
  