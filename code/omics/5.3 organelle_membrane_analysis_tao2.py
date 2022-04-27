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




# analyze all membranes as a whole
# bar plot
value_df = membrane_only_ratio
column_select1 = column_select10


g1 = value_df[value_df["sample_ID"]==5]
g1 = g1[column_select1]
g10 = g1.mean()
g100 = pd.DataFrame(g10)
g100.columns = ["C:N=5"]

g2 = value_df[value_df["sample_ID"]==115]
g2 = g2[column_select1]
g20 = g2.mean()
g200 = pd.DataFrame(g20)
g200.columns = ["C:N=115"]

# combine
pd_null = pd.concat([g100, g200], axis=1)
pd_null["compartment"] = list(pd_null.index)
pd_null["Relative change"] = 2*(pd_null["C:N=115"]-pd_null["C:N=5"])/(pd_null["C:N=5"]+pd_null["C:N=115"])*100
#pd_null["fold_change"] = pd_null["C:N=50"] / pd_null["C:N=5"]
pd_null = pd_null.sort_values(by=['Relative change'], ascending=False)

# bar plot
x0='compartment'
y0='Relative change'
plt.figure()
sns.set_style('darkgrid')
sns.barplot(x=x0, y=y0, data=pd_null, capsize=.2)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0 + " in C:N=115 vs C:N=5 (%)", fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/membrane area ratio correlation analysis2.pdf", bbox_inches='tight')



