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








# make a standard graph 4.94*4.54 device size
fit1 <- lm(WT_4 ~ WT_0, data = df)  
ggplot(df, aes(x=WT_0 , y=WT_4, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 4h") +
  xlim(0, 0.45) + ylim(0,0.45) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.35), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3),
                           " \nP value =",signif(summary(fit1)$coef[2,4], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(WT_12 > 0.05, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)





fit1 <- lm(WT_8 ~ WT_0, data = df)  
ggplot(df, aes(x=WT_0 , y=WT_8, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 8h") +
  xlim(0, 0.45) + ylim(0,0.45) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.35), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3),
                           " \nP value =",signif(summary(fit1)$coef[2,4], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(WT_12 > 0.05, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)




fit1 <- lm(WT_12 ~ WT_0, data = df)  
ggplot(df, aes(x=WT_0 , y=WT_12, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Zinc limitation 0h",
       y = "Zinc limitation 12h") +
  xlim(0, 0.45) + ylim(0,0.45) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.35), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3),
                           " \nP value =",signif(summary(fit1)$coef[2,4], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(WT_12 > 0.05, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)























df$WT_8_vs_WT_0 <- (df$WT_8-df$WT_0)/df$WT_0
df$WT_8_vs_WT_0_abs <- abs(df$WT_8_vs_WT_0)

df_filter_mass <- df[df$WT_0 >=0.01, ]
df_filter_mass <- df_filter_mass[df_filter_mass$WT_8_vs_WT_0_abs >=0.2, ] 

ggplot(df_filter_mass, aes(x=reorder(compartment, WT_8_vs_WT_0), y=WT_8_vs_WT_0)) + 
  geom_bar(position="dodge", stat="identity", fill = "red", alpha=0.4) +
  xlab("Cellular Component") + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12))+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position="none")+
  coord_flip()


