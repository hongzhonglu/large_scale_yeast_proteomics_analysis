# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns

# input the datasets
data_rosemary = pd.read_excel("data/proteomics/organell_protein_ratio_rosemary.xlsx")
data_jianye = pd.read_excel("data/proteomics/organell_protein_ratio_jianye.xlsx")
data_tao2 = pd.read_excel("data/proteomics/organell_protein_ratio_tao2.xlsx")


# set the unified column name
colname0 = list(data_jianye.columns)[1:]
colname0 = [x for x in colname0 if x !="total_pro_volume"]


# unify other datasets
data_jianye0 = data_jianye[colname0]
s1 = data_tao2["sample_ID"].tolist()
data_tao2["sample_ID"] = ["C_N_" + str(x) for x in s1]
data_tao2 = data_tao2[colname0]
# calculate the average values

s2 = data_rosemary["dilution rate (/h)"].tolist()
data_rosemary["sample_ID"] = ["D=" + str(x) for x in s2]
data_rosemary = data_rosemary[colname0]


# calculate the average values
data_tao20 = data_tao2.groupby(['sample_ID']).mean()
data_tao20["sample_ID"] = list(data_tao20.index)

data_rosemary0 = data_rosemary.groupby(['sample_ID']).mean()
data_rosemary0["sample_ID"] = list(data_rosemary0.index)


# calculate the relative change
data_jianye_t = data_jianye0.transpose()
data_jianye_t.columns = data_jianye_t.iloc[139,:]
data_jianye_t = data_jianye_t.iloc[0:139,6:9]

data_rosemary_t =data_rosemary0.transpose()
data_rosemary_t.columns = data_rosemary_t.iloc[139,:]
data_rosemary_t = data_rosemary_t.iloc[0:139,3:6]

data_tao2_t = data_tao20.transpose()
data_tao2_t.columns = data_tao2_t.iloc[139,:]
data_tao2_t = data_tao2_t.iloc[0:139,:]
col_order = ['C_N_5','C_N_30', 'C_N_50', 'C_N_115']
data_tao2_t = data_tao2_t[col_order]

# next
def calculate_relative_change(data_t):
    data_t0 = data_t.copy()
    all_col = list(data_t0.columns)
    ref = all_col[0]
    for x in all_col:
        print(x)
        data_t0[x] = 2*(data_t[x] - data_t[ref]) / (data_t[ref]+data_t[x])
    return data_t0.iloc[:,1:]

data_jianye_t0 = calculate_relative_change(data_jianye_t)
data_rosemary_t0 = calculate_relative_change(data_rosemary_t)
data_tao2_t0 = calculate_relative_change(data_tao2_t)

# combine
pd_combine = pd.concat([data_jianye_t0, data_rosemary_t0, data_tao2_t0], axis=1)
pd_combine = pd_combine.dropna()
pd_combine.to_excel("data/proteomics/organelle_ratio_combine.xlsx")



"""
# pheatmap in R
library(pheatmap)
library(readxl)
organelle_ratio_combine <- read_excel("organelle_ratio_combine.xlsx")
rownames(organelle_ratio_combine) <- organelle_ratio_combine$organelle
organelle_ratio_combine0 <-  subset(organelle_ratio_combine, select = -c(organelle, C_N_30))
rownames(organelle_ratio_combine0) <- organelle_ratio_combine$organelle
dat <- cbind(matrix(rnorm(120), 30, 40), matrix(sample(15, 120, T), 30))
my.breaks <- c(seq(-1.5, 0, by=0.1), seq(0.1, 1.5, by=0.1)) 
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "orange", "red", "purple"))(length(my.breaks)/2))
pheatmap(organelle_ratio_combine0,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = TRUE,
         cluster_col = TRUE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 3,
         color = my.colors,
         breaks = my.breaks)
"""
