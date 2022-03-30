# gene level analysis

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

# input the dn_ds value of sce
dn_ds_sce = pd.read_csv("data/dn_ds_sce.csv")


# just select one column
data_carl = protein_copy_all1[['gene', 'mmol/gDW_carl']]
data_carl = data_carl.dropna()
df_carl = pd.merge(left=data_carl, right=pro_size, left_on=['gene'], right_on=['locus'], how="left")
df_carl = df_carl.sort_values(by=['Total_Volume'], ascending=False)
df_carl["dN_dS"] = singleMapping(dn_ds_sce["dN_dS"],dn_ds_sce["locus"],df_carl["gene"])

# plot
x0='Total_Volume'
y0='mmol/gDW_carl'
fig, ax = plt.subplots(figsize=(6, 6))
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=df_carl)

sns.kdeplot(
    data=df_carl,
    x=x0,
    y=y0,
    levels=5,
    fill=True,
    alpha=0.6,
    cut=2,
    ax=ax,
)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# select two datasets
data_two_sample = protein_copy_all1[['gene','Glucose_phase(mmol/gDW)','mmol/gDW_carl']]
data_two_sample = data_two_sample.dropna()

data_two_sample["log1"] = np.log10(data_two_sample['Glucose_phase(mmol/gDW)'])
data_two_sample["log2"] = np.log10(data_two_sample['mmol/gDW_carl'])

# plot
x0='log1'
y0='log2'
plt.figure()
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=data_two_sample)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.plot([0, 8], [0, 8], linewidth=2)

# add a reference for the batch conditions
data_three_sample = protein_copy_all1[['gene','Glucose_phase(mmol/gDW)','mmol/gDW_carl', 'Min_ave_aerobic(mmol/gDW)']]
df = data_three_sample.describe()





