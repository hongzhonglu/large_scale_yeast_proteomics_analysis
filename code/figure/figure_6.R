library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)
library(plotly)
library(Rtsne) # tSNE in an acronym for t-Distributed Neighbor Embedding is a statistical method that is mainly used to visualize high-dimensional data
library(umap) # umap is similar to tSNE, but more efficient
# part 1 proteomics analysis-main part
# data and sample
combine <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/mass_fraction_NCB_for_yeast_IO.xlsx")


# analyze the intersection of all samples
na_counts_per_row <- rowSums(is.na(combine[,2:27]))
# If you want the result as a dataframe
na_counts_df <- data.frame(row_NA_count = na_counts_per_row, row.names = combine$gene)
gene_remove <- rownames(na_counts_df)[which(na_counts_df$row_NA_count >=20)]

combine <- combine[!(combine$gene %in% gene_remove),]

physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")
physiology_collection0 <- physiology_collection[!duplicated(physiology_collection$condition_unique),]
combine0 <- combine[, colnames(combine) %in% physiology_collection0$sampleID]

# input the compartment annotation
compartment_annotation <- read_excel("/Users/xluhon/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/nature_chemical_biology_datatset_2024/IO_gene_compartment_filter.xlsx")
#compartment_annotation['fold'] <- compartment_annotation$annotation_curation / compartment_annotation$annotation_combine
# remove 'membrane'? as it actually covers "membrane protein" from different organelles
compartment_annotation <- compartment_annotation[compartment_annotation$compartment !="membrane",]
compartment0 <- compartment_annotation[compartment_annotation$compartment !="mitochondrion_unassigned",]
compartment0 <- compartment0[compartment0$compartment !="cytoplasm",]
compartment_sum <- as.data.frame(table(compartment0$compartment))
compartment_sum <- compartment_sum[order(compartment_sum[[2]], decreasing = TRUE), ]
source_sum <- as.data.frame(table(compartment0$source))

# bar plot to analyze the protein num in each main component
ggplot(data=compartment_sum[1:10,], aes(x = reorder(Var1, Freq), y = Freq)) + 
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


ggplot(data=source_sum, aes(x = reorder(Var1, Freq), y = Freq)) + 
  geom_bar(stat="identity",fill="#619CFF") +
  xlab("Annotation source") + 
  ylab("Protein count") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))




# calculate the organelle mass fraction variance 
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_NCB_yeast_IO.xlsx")
ProMassRatio <- ProMassRatio[,2:27]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
sd0 <- apply(subset(ProMassRatio, select = 2:26), 1, sd, na.rm=TRUE) 
mean0 <- apply(subset(ProMassRatio, select = 2:26), 1, mean, na.rm=TRUE) 
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
ProMassRatio_ss <- ProMassRatio1[, colnames(ProMassRatio1) %in% physiology_collection0$sampleID] # remove the duplicated ones
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
  xlim(0.95, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "Correlation coefficient between samples",
       y = "Density") 



# Heatmap and PCA plot for IO under all conditions
# heatmap
library("pheatmap")
ProMassRatio_IO <- ProMassRatio1[,c(2:26)]
ProMassRatio_IO100 <- as.matrix(ProMassRatio_IO)
#heatmap(ProMassRatio_IO)
pheatmap(ProMassRatio_IO100, scale="row",
         show_colnames = TRUE,
         show_rownames = FALSE,
         border_color=NA)
# PCA plot
combine11 <- t(ProMassRatio_IO100)
combine11[is.na(combine11)] <- 0 # here NA value was replaced as 0

iris.umap = umap(combine11, n_components = 2, random_state = 15) 

layout <- iris.umap[["layout"]] 
layout <- data.frame(layout) 
label11 <- rownames(layout)
final <- cbind(layout, label11) 

fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~label11, type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Source')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig 


ProMassRatio_IO$compartment <- ProMassRatio1$compartment
# check the resource allocation under different conditions: with oxygen, without oxygen, with inhibition
ProMassRatio_IO2 <- ProMassRatio_IO[,str_detect(colnames(ProMassRatio_IO), "O2")|str_detect(colnames(ProMassRatio_IO), "antimycin")]
ProMassRatio_IO2$compartment <- ProMassRatio_IO$compartment
#combine the same condition
condition <- str_trim(colnames(ProMassRatio_IO2), side = "both")
condition <- str_replace_all(condition, "_[:digit:]", "")

colnames(ProMassRatio_IO2) <- condition
ProMassRatio_IO2 <- ProMassRatio_IO2[!is.na(ProMassRatio_IO2$`IO_SD108_glucose + O2`),]
Pro_mass_select0 <- ProMassRatio_IO2 [, !colnames(ProMassRatio_IO2 ) %in% c("compartment")]

# calculate the average of columns
df <- as.data.frame(sapply(split.default(Pro_mass_select0, names(Pro_mass_select0)), rowMeans))
df$compartment <- ProMassRatio_IO2$compartment
df <- df[df$compartment !="cytoplasm", ]
df <- df[df$compartment !="mitochondrion_unassigned", ]

colnames(df) <- c("IO_no_O2",       "IO_antimycin",     "IO_O2",  "compartment" ) 

fit1 <- lm( IO_no_O2 ~ IO_O2, data = df)  
ggplot(df, aes(x=IO_O2 , y=IO_no_O2, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "IO_O2",
       y = "IO_no_O2 ") +
  xlim(0, 0.27) + ylim(0,0.27) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.2), hjust = 0, 
             label = paste("Adj R² = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(IO_no_O2 > 0.1, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)



fit1 <- lm( IO_antimycin ~ IO_O2, data = df)  
ggplot(df, aes(x=IO_O2 , y=IO_antimycin, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "IO_O2",
       y = "IO_antimycin ") +
  xlim(0, 0.27) + ylim(0,0.27) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.2), hjust = 0, 
             label = paste("Adj R² = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(IO_antimycin > 0.1, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)






#############################################################################################3
# Heatmap and PCA plot for sce
ProMassRatio_sce_all <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio_sce1 <- ProMassRatio_sce_all[, str_detect(colnames(ProMassRatio_sce_all), 'sce_FY4')|str_detect(colnames(ProMassRatio_sce_all), 'sce_CEN.PK')]
ProMassRatio_sce100 <- as.matrix(ProMassRatio_sce1)
#heatmap(ProMassRatio_sce1)
pheatmap(ProMassRatio_sce100, scale="row",
         show_colnames = TRUE,
         show_rownames = FALSE,
         border_color=NA)

# PCA plot
combine11 <- t(ProMassRatio_sce100)
combine11[is.na(combine11)] <- 0 # here NA value was replaced as 0

iris.umap = umap(combine11, n_components = 2, random_state = 15) 

layout <- iris.umap[["layout"]] 
layout <- data.frame(layout) 
label11 <- rownames(layout)
final <- cbind(layout, label11) 

fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~label11, type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='Source')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 
fig 





# check the resource allocation under different conditions: with oxygen, without oxygen, with inhibition
ProMassRatio_sce2 <- ProMassRatio_sce1[,str_detect(colnames(ProMassRatio_sce1), "O2")|str_detect(colnames(ProMassRatio_sce1), "antimycin")]
ProMassRatio_sce2$compartment <- ProMassRatio_sce_all$compartment
#combine the same condition
condition <- str_trim(colnames(ProMassRatio_sce2), side = "both")
condition <- str_replace_all(condition, "_[:digit:]", "")

colnames(ProMassRatio_sce2) <- condition
ProMassRatio_sce2 <- ProMassRatio_sce2[!is.na(ProMassRatio_sce2$`sce_CEN.PK_glucose + O2`),]
Pro_mass_select0 <- ProMassRatio_sce2 [, !colnames(ProMassRatio_sce2 ) %in% c("compartment")]

# calculate the average of columns
df <- as.data.frame(sapply(split.default(Pro_mass_select0, names(Pro_mass_select0)), rowMeans))
df$compartment <- ProMassRatio_sce2$compartment
df <- df[df$compartment !="cytoplasm", ]
df <- df[df$compartment !="mitochondrion_unassigned", ]

colnames(df) <- c("sce_no_O2",       "sce_antimycin",     "sce_O2",  "compartment" ) 

fit1 <- lm( sce_no_O2 ~ sce_O2, data = df)  
ggplot(df, aes(x=sce_O2 , y=sce_no_O2, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "sce_O2",
       y = "sce_no_O2 ") +
  xlim(0, 0.45) + ylim(0,0.45) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.35), hjust = 0, 
             label = paste("Adj R² = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(sce_no_O2 > 0.1, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)



fit1 <- lm( sce_antimycin ~ sce_O2, data = df)  
ggplot(df, aes(x=sce_O2 , y=sce_antimycin, label=compartment)) +
  geom_point(size=4, shape=1,colour='#E69F00') +
  geom_smooth(method=lm) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  labs(x = "sce_O2",
       y = "sce_antimycin ") +
  xlim(0, 0.45) + ylim(0,0.45) +
  geom_abline(slope=1, intercept=0, linetype=2, size=1.5, colour = "grey") +
  geom_label(aes(x = 0, y = 0.35), hjust = 0, 
             label = paste("Adj R² = ",signif(summary(fit1)$adj.r.squared, 3),
                           "\nIntercept =",signif(fit1$coef[[1]],3),
                           " \nSlope =",signif(fit1$coef[[2]], 3)),
             label.size = NA)+
  geom_text(aes(label=ifelse(sce_antimycin > 0.1, as.character(compartment),'')),hjust=-0.1,vjust=-0.1)





# Supplementary file
# stacked barplot
#organelle_s <- c('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum', 'fungal-type vacuole','peroxisome','ribosome')
ProMassRatio2 <- ProMassRatio1[ProMassRatio1$compartment %in% organelle, ]
test <- ProMassRatio2[, str_detect(colnames(ProMassRatio2), 'IO_SD108_C')]
test$compartment <- ProMassRatio2$compartment
long_DF <- test %>% gather(growth, mass_fraction, 1:5)

long_DF %>% 
  group_by(growth, compartment) %>% 
  ggplot(aes(x = growth, y = mass_fraction, group = compartment, fill = compartment)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# plus sce
# calculate the organelle mass fraction variance 
ProMassRatio_sce_all <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio_sce0 <- ProMassRatio_sce_all[, str_detect(colnames(ProMassRatio_sce_all), 'sce_FY4_C')]
ProMassRatio_sce0$compartment <- ProMassRatio_sce_all$compartment
ProMassRatio_sce01 <- ProMassRatio_sce0[ProMassRatio_sce0$compartment %in% organelle, ]

long_DF <- ProMassRatio_sce01 %>% gather(growth, mass_fraction, 1:4)

long_DF %>% 
  group_by(growth, compartment) %>% 
  ggplot(aes(x = growth, y = mass_fraction, group = compartment, fill = compartment)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))












