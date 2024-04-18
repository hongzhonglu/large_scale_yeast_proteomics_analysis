library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)
# corrplot
# part 1 proteomics analysis-main part
# data and sample
library(readxl)
mass_fraction_combine <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/mass_fraction_combine.xlsx")
physiology_collection <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/physiology_collection.xlsx")
physiology_collection0 <- physiology_collection[!duplicated(physiology_collection$condition_unique),]

mass_fraction_combine1 <- mass_fraction_combine[, colnames(mass_fraction_combine) %in% physiology_collection0$sampleID]



combine0 <- mass_fraction_combine1

M <- cor(combine0, method = "pearson", use = "pairwise.complete.obs")
corrplot(M[1:50, 1:50], tl.pos='n', type = "upper", shade.lwd=0.0001)

ss <- as.vector(M[upper.tri(M)])
df <- data.frame(value=ss)

# plot density plot
ggplot(df, aes(value)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "black") +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(alpha = 0.5, fill = "lightgreen")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  labs(x = "R2",
       y = "Density") 






# part2 protein compartment







# part3 protein properties - mass, volume, sectional area








