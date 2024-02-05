# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
# input the membrane size data
# volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation.xlsx") # not curated
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate


volume_size_tr = volume_size.transpose()
volume_size_tr = volume_size_tr.iloc[1:, :]
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[0])



column_select1 = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet',
                  'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'cytosolic ribosome','mitochondrial ribosome','nucleolus']
column_select1 = [x for x in column_select1 if "ribosome" not in x]

volume_size_tr0 = volume_size_tr0.iloc[1:, :]
volume_size_tr0 = volume_size_tr0[column_select1]


# enzyme
enzyme_size = pd.read_excel("data/proteomics/GEM_volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate
enzyme_size_tr = enzyme_size.transpose()

enzyme_size_tr0 = enzyme_size_tr.rename(columns=enzyme_size_tr.iloc[0])
enzyme_size_tr0 = enzyme_size_tr0.iloc[1:, :]
enzyme_size_tr0 = enzyme_size_tr0[column_select1]

# calculate the ratio
ratio_df = enzyme_size_tr0.copy()
for x in column_select1:
    ratio_df[x] = enzyme_size_tr0[x]/volume_size_tr0[x]

ratio_df["id"] = list(ratio_df.index)

ratio_df["dilution rate (/h)"] = singleMapping(physiology_data["dilution rate (/h)"], physiology_data["kinetic"],ratio_df["id"])
ratio_df = ratio_df[column_select1+["dilution rate (/h)"]]



s2 = ratio_df["dilution rate (/h)"].tolist()
colname0 = list(ratio_df.columns)
colname0 = [x for x in colname0 if x !="dilution rate (/h)"]
# change the data format
ratio_df0 = ratio_df[colname0]
for x in colname0:
    ratio_df0[x] = pd.to_numeric(ratio_df0[x])

ratio_df0["sample_ID"] = ["D=" + str(x) for x in s2]

ratio_df01 = ratio_df0.groupby(['sample_ID']).mean()

ratio_df01["dilution rate (/h)"] = list(dict.fromkeys(s2))



# plot the figure
sns.set_style("darkgrid")
plt.figure()
sns.lineplot(x='dilution rate (/h)', y='value', hue='variable', style="variable",
             data=pd.melt(ratio_df01, ['dilution rate (/h)']))
#plt.axvline(x=0.18, color='g', linestyle='--')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('result/figure/rose_miu_metabolic_enzyme_per_organelle.pdf', bbox_inches='tight')