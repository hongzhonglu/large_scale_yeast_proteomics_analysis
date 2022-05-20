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

# calculate the ratio
data_tao2["c_to_m"] = data_tao2["cytosol"]/data_tao2["mitochondrion"]
data_jianye0["c_to_m"] = data_jianye0["cytosol"]/data_jianye0["mitochondrion"]
data_rosemary["c_to_m"] = data_rosemary["cytosol"]/data_rosemary["mitochondrion"]


data_rosemary = data_rosemary[["dilution rate (/h)", "c_to_m"]]
data_jianye0 = data_jianye0[["sample_ID", "c_to_m"]]
data_jianye0["sample_ID"] = data_jianye0["sample_ID"].str.replace("D=","").str.replace("_M","")
data_jianye0.sample_ID = pd.to_numeric(data_jianye0.sample_ID, errors='coerce')


x0="dilution rate (/h)"
y0="c_to_m"
title0 = 'result/figure/rosemary_ratio_' + y0 + '.pdf'
plt.figure(figsize=(4,4))
sns.lineplot(x=x0, y=y0, data=data_rosemary, marker="o")
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0, 0.4)
plt.axvline(x=0.18, color='k', linestyle='--')
plt.savefig(title0, bbox_inches='tight')


x0="sample_ID"
y0="c_to_m"
title0 = 'result/figure/jianye_ratio_' + y0 + '.pdf'
plt.figure(figsize=(4,4))
sns.lineplot(x=x0, y=y0, data=data_jianye0, marker="o")
plt.xlabel("dilution rate (/h)", fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0, 0.4)
plt.axvline(x=0.284, color='k', linestyle='--')
plt.savefig(title0, bbox_inches='tight')









data_tao2 = data_tao2[["sample_ID", "c_to_m"]]
data_tao2.sample_ID = data_tao2.sample_ID.str.replace("C_N_", "")
data_tao2.sample_ID = pd.to_numeric(data_tao2.sample_ID, errors='coerce')

x0="sample_ID"
y0="c_to_m"
title0 = 'result/figure/tao2_ratio_' + y0 + '.pdf'
print(title0)
plt.figure(figsize=(4,4))
sns.barplot(x=x0, y=y0, data=data_tao2, capsize=.2)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.savefig(title0, bbox_inches='tight')



