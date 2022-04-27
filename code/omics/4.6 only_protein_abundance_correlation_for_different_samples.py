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

# select two datasets
data_two_sample = protein_copy_all1[['gene','Glucose_phase(mmol/gDW)','mmol/gDW_carl']]
data_two_sample = data_two_sample.dropna()

data_two_sample["abs1"] = (data_two_sample['Glucose_phase(mmol/gDW)'])
data_two_sample["abs2"] = (data_two_sample['mmol/gDW_carl'])

from scipy.stats import pearsonr
corr1 = pearsonr(data_two_sample["abs1"], data_two_sample["abs2"])




# plot
x0='abs1'
y0='abs2'
plt.figure()
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=data_two_sample)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# add a reference for the batch conditions
data_three_sample = protein_copy_all1[['gene','Glucose_phase(mmol/gDW)','mmol/gDW_carl', 'Min_ave_aerobic(mmol/gDW)']]
df = data_three_sample.describe()


# how the result is going when considering the volume of proteins?
data_two_sample["protein_volume"] = singleMapping(pro_size["Total_Volume"], pro_size["locus"],data_two_sample["gene"])
data_two_sample["v1"] = (data_two_sample["Glucose_phase(mmol/gDW)"] * data_two_sample["protein_volume"])
data_two_sample["v2"] = (data_two_sample["mmol/gDW_carl"] * data_two_sample["protein_volume"])
# plot
x0='v1'
y0='v2'
plt.figure()
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=data_two_sample)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
from scipy.stats import pearsonr
corr2, _ = pearsonr(data_two_sample["v1"], data_two_sample["v2"])