# gene level analysis

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the pro structure size data
pro_size = pd.read_excel("data/other_species/ecoli_structure_size_part1.xlsx")

pro_size = pro_size[['Protein','Total_Volume']]
pro_size["Protein"] = pro_size["Protein"].str.replace("-F1-model_v2","").str.replace("AF-","")

# for all the proteins
sns.displot(pro_size, x="Total_Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["Total_Volume"]]
df1 = df.describe()

# further input the protein length and molecular weight
pro_info = pd.read_csv("data/other_species/ecoli.txt", sep="\t")
pro_size["MW"] = singleMapping(pro_info["Mass"], pro_info["Entry"],pro_size["Protein"])
pro_size["pro_length"] = singleMapping(pro_info["Length"], pro_info["Entry"],pro_size["Protein"])

pro_size["MW"] = pro_size["MW"].str.replace(",", "")
pro_size.MW = pd.to_numeric(pro_size.MW, errors='coerce')
a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
# filter by protein length to remove too short proteins

pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]
pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)

# one interesting idea is used sce formula to calculate the protein volume in other species
pro_size["calculated_volume"] = 1.06019171e-03*pro_size["MW"] - 1.10587455
# compare the predicted and calculated for e.coli
x0 = "MW"
y0 = "Total_Volume"
y1 = "calculated_volume"
plt.figure(figsize=(4, 4))
sns.lineplot(x=x0, y=y1, data=pro_size)
sns.scatterplot(x=x0, y=y0, data=pro_size)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y1, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)








# spo
# input the pro structure size data
pro_size = pd.read_excel("data/other_species/spo_structure_size_part1.xlsx")

pro_size = pro_size[['Protein','Total_Volume']]
pro_size["Protein"] = pro_size["Protein"].str.replace("-F1-model_v2","").str.replace("AF-","")

# for all the proteins
sns.displot(pro_size, x="Total_Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["Total_Volume"]]
df1 = df.describe()

# further input the protein length and molecular weight
pro_info = pd.read_csv("data/other_species/spo.txt", sep="\t")
pro_size["MW"] = singleMapping(pro_info["Mass"], pro_info["Entry"],pro_size["Protein"])
pro_size["pro_length"] = singleMapping(pro_info["Length"], pro_info["Entry"],pro_size["Protein"])

pro_size["MW"] = pro_size["MW"].str.replace(",", "")
pro_size.MW = pd.to_numeric(pro_size.MW, errors='coerce')
a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
# filter by protein length to remove too short proteins

pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]
pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)
sns.displot(pro_size, x="volume_per_kda", stat="density", common_norm=False)


# one interesting idea is used sce formula to calculate the protein volume in other species
pro_size["calculated_volume"] = 1.06019171e-03*pro_size["MW"] - 1.10587455
# compare the predicted and calculated for e.coli
x0 = "MW"
y0 = "Total_Volume"
y1 = "calculated_volume"
plt.figure(figsize=(4, 4))
sns.lineplot(x=x0, y=y1, data=pro_size)
sns.scatterplot(x=x0, y=y0, data=pro_size)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y1, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)









# import Transcription factor
from scipy.stats import ttest_ind
pro_size00 = pro_size.copy()
TRN_spo = pd.read_csv("data/other_species/spo_TF_list.tsv", sep="\t")
pro_size00["TF"] = None
pro_size00["TF"][pro_size00["Protein"].isin(TRN_spo["UniProt ID"])] = "Yes"
pro_size00["TF"][~pro_size00["Protein"].isin(TRN_spo["UniProt ID"])] = "No"
pro_g1 = pro_size00[pro_size00["TF"]=="Yes"]
pro_g2 = pro_size00[pro_size00["TF"]=="No"]
pro_g2 = pro_g2[pro_g2["pro_length"] >= min(pro_g1["pro_length"])]
pro_g2 = pro_g2[pro_g2["pro_length"] <= max(pro_g1["pro_length"])]
# combine two pandas
pro_c = pd.concat([pro_g1, pro_g2], axis=0)
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])
