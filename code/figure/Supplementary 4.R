library(readr)
library(readxl)
library(tidyverse)
library(hongR)
library(corrplot)
library(ggplot2)
library(dplyr)

# compartment
ProMassRatio <- read_excel("~/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/proteomics/ProMassRatio_across_compartment_combine.xlsx")
ProMassRatio <- ProMassRatio[,2:277]
ProMassRatio[ProMassRatio <0.0000000000001] <- NA
ProMassRatio1 <- ProMassRatio[ProMassRatio$compartment !="cytoplasm", ]
ProMassRatio1 <- ProMassRatio1[ProMassRatio1$compartment !="mitochondrion_unassigned", ]





# get the protein abundance under batch condition

batch_glucose <- ProMassRatio1 %>% select(1:4)
batch_glucose$row_mean <- rowMeans(batch_glucose[, 2:4])


batch_ethanol <- ProMassRatio1 %>% select(8:10)
batch_ethanol$row_mean <- rowMeans(batch_ethanol[, 1:3])

# combine

df <- batch_glucose[,c(1,5)]
colnames(df) <- c("compartment","batch_glucose")
df$batch_ethanol <- batch_ethanol$row_mean
df <- df[!is.na(df$batch_glucose),]
#large components with PMF greater than 0.0035; medium components with PMF between 0.00035-0.0035; small components with PMF smaller than 0.00035 

df1 <- df[df$batch_glucose >= 0.0035,]

df2 <- df[df$batch_glucose < 0.0035 & df$batch_glucose >=0.00035,]

df3 <- df[df$batch_glucose < 0.00035,]


# scatter plot
ggplot(df1, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0.0035, 0.4) +
  ylim(0.0035, 0.4) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df1)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6




# scatter plot
ggplot(df2, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0.00035, 0.0035) +
  ylim(0.00035, 0.0035) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df2)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6






# scatter plot
ggplot(df3, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0, 0.00035) +
  ylim(0, 0.00035) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df3)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6


















#######
# Chemostat condition
#######

batch_glucose <- ProMassRatio1 %>% select(c("prot.1","prot.2", "prot.3"))
batch_glucose$row_mean <- rowMeans(batch_glucose[, 1:3])


batch_ethanol <- ProMassRatio1 %>% select(c("prot.19","prot.20", "prot.21"))
batch_ethanol$row_mean <- rowMeans(batch_ethanol[, 1:3])



# combine

df <- batch_glucose[,c(1,4)]
colnames(df) <- c("compartment","batch_glucose")
df$batch_ethanol <- batch_ethanol$row_mean
df <- df[!is.na(df$batch_glucose),]
#large components with PMF greater than 0.0035; medium components with PMF between 0.00035-0.0035; small components with PMF smaller than 0.00035 

df1 <- df[df$batch_glucose >= 0.0035,]

df2 <- df[df$batch_glucose < 0.0035 & df$batch_glucose >=0.00035,]

df3 <- df[df$batch_glucose < 0.00035,]


# scatter plot
ggplot(df1, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0.0035, 0.3) +
  ylim(0.0035, 0.3) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df1)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6




# scatter plot
ggplot(df2, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0.00035, 0.0035) +
  ylim(0.00035, 0.0035) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df2)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6






# scatter plot
ggplot(df3, aes(x=batch_glucose, y=batch_ethanol)) +
  geom_point(colour = "black", size = 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  labs(x = "Mass fraction of components",
       y = "Mass fraction of components") +
  xlim(0, 0.00035) +
  ylim(0, 0.00035) +
  geom_smooth(method=lm , color="red", se=FALSE)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))

# 拟合线性模型
model <- lm(batch_ethanol~ batch_glucose, data = df3)
# 提取R²值
r_squared <- summary(model)$r.squared
print(paste("R² =", round(r_squared, 4)))  # 输出：R² = 0.6







