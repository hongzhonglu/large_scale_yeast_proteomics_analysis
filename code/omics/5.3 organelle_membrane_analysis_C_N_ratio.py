# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *
import seaborn as sns


# input the membrane size data
# membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation.xlsx") # not curated
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_tao2.xlsx") # curated based on cell size under different growth rate

membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])




# only take rosemery physiology dataset

# only take rosemery proteomics
membrane_size_1 = membrane_size_tr0[membrane_size_tr0.index.str.contains("fmol")]
membrane_size_1["sample_ID"] = [5,5,115,115,30,30,50,50]
# combine the physiological datasets and proteomics datasets
combine_data = membrane_size_1


# check the ratio of each organelle membrane size
column_select = list(combine_data.columns)
column_select1 = [x for x in column_select if "membrane" in x]
column_select10 = [x for x in column_select1 if "component" not in x]
column_select10 = [x for x in column_select10 if "contact site" not in x]
column_select10 = [x for x in column_select10 if "raft" not in x]
column_select10 = [x for x in column_select10 if x!="membrane"]
column_select10 = [x for x in column_select10 if "network" not in x]
column_select10 = [x for x in column_select10 if "space" not in x]

# because the size of plasma membrane and fungal-type vacuole membrane is so bigger, exclude them firstly
# column_select10 = [x for x in column_select10 if x!="plasma membrane"]
# column_select10 = [x for x in column_select10 if x!="fungal-type vacuole membrane"]
column_select10 = [x for x in column_select10 if x!="mitochondrial membrane"] # remove duplicates?
column_select10 = [x for x in column_select10 if x!="vacuolar membrane"] # remove duplicates?


membrane_only = combine_data[column_select10]

# calculate the ratio of each organelle membrane relative to total membrane
membrane_only_ratio = membrane_only.copy()
i0 = -1
for i, x in membrane_only.iterrows():
    print(i)
    i0=i0+1
    ss = list(x)
    for j in range(len(ss)):
        print(j)
        ratio = ss[j]/sum(ss)
        membrane_only_ratio.iloc[i0,j] = ratio

membrane_only_ratio["sample_ID"] = combine_data["sample_ID"]


# plot the figure in occupied area
x0 = 'sample_ID'
for y0 in column_select10:
    title0 = 'result/figure/tao2_' + y0 + '.pdf'
    print(title0)
    plt.figure()
    sns.barplot(x=x0, y=y0, data=combine_data, capsize=.2)
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied area",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(rotation=90)
    plt.savefig(title0, bbox_inches='tight')


# plot the figure in ratio
x0 = "sample_ID"
for y0 in column_select10:
    title0 = 'result/figure/tao2_ratio_' + y0 + '.pdf'
    print(title0)
    plt.figure()
    sns.barplot(x=x0, y=y0, data=membrane_only_ratio, capsize=.2)
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied ratio",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(rotation=90)
    plt.savefig(title0, bbox_inches='tight')