# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Kate
sample_info_kate = pd.read_excel("data/proteomics/datasets_kate_2020.xlsx", sheet_name="sample_information")
sample_info_kate = sample_info_kate.transpose()
sample_info_kate1 = sample_info_kate.iloc[1: , :]
sample_info_kate1.columns = sample_info_kate.iloc[0]
sample_info_kate1["sample_ID"] = list(sample_info_kate1.index)


# input the membrane size data
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_kate.xlsx") # curated based on cell size under different growth rate


volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])
volume_size_tr0 = volume_size_tr0.iloc[2: ,:]


# check other volume
column_select = list(volume_size_tr0.columns)
column_select1 = [x for x in column_select if "membrane" not in x]
column_select1 = [x for x in column_select1 if "wall" not in x]
column_select1 = [x for x in column_select1 if "site" not in x]
column_select1 = [x for x in column_select1 if "tip" not in x]
column_select1 = [x for x in column_select1 if "pore" not in x]
column_select1 = [x for x in column_select1 if "catalytic" not in x]

column_select1 = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet',
                  'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'cytosolic ribosome','mitochondrial ribosome','nucleolus']

# note: The nucleolus is a region found within the cell nucleus that is concerned with producing and assembling the cell's ribosomes.
# 'mitochondrial ribosome'
combine_data2 = volume_size_tr0[column_select1]
combine_data2["sample_ID"] = sample_info_kate1["sample_ID"]
combine_data2['Assigned cell cycle phase '] = sample_info_kate1['Assigned cell cycle phase ']
combine_data2["sample_ID"] = combine_data2['Assigned cell cycle phase '] + "@" + combine_data2["sample_ID"]


# calculate the average value
all_ID = combine_data2["sample_ID"].to_list()
g1 = all_ID[0:18]
g1_n = [x.split("@")[0] + "_1" for x in g1]
g2 = all_ID[18:27]
g2_n = [x.split("@")[0] + "_2" for x in g2]
g3 = all_ID[27:]
g3_n = [x.split("@")[0] + "_3" for x in g3]
all_ID_n = g1_n + g2_n + g3_n
combine_data2["sample_ID"] = all_ID_n


# plot the figure in ratio
x0 = "sample_ID"
for y0 in column_select1:
    title0 = 'result/figure/kate_volume_' + y0 + '.pdf'
    print(title0)
    plt.figure()
    sns.barplot(x=x0, y=y0, data=combine_data2, capsize=.2)
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " volume",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(rotation=90)
    plt.savefig(title0, bbox_inches='tight')