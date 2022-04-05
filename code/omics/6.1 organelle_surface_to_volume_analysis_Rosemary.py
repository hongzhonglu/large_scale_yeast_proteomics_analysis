# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]

# input the volume size data
# volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation.xlsx") # not curated
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate
volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])
volume_size_tr0 = volume_size_tr0.iloc[2:20,:]

column_select1 = ['mitochondrion', 'nucleus', 'endoplasmic reticulum','endosome',
                  'fungal-type vacuole','peroxisome','Golgi apparatus']
#volume_size_tr0 = volume_size_tr0[column_select1]
#volume_size_tr0["sample_ID"] = list(volume_size_tr0.index)


# input the membrane size data
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate
membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])
membrane_size_tr0 = membrane_size_tr0.iloc[2:20,:]



# check the ratio of each organelle membrane size
column_select = list(membrane_size_tr0.columns)
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
membrane_not_use = ['prospore membrane', 'cellular bud membrane','late endosome membrane', 'plasma membrane']
column_select10 = [x for x in column_select10 if x not in membrane_not_use]


membrane_size_tr0 = membrane_size_tr0[column_select10]


# example analysis
organelle1 = ['mitochondrion']
organelle1_m = ['mitochondrial outer membrane',  'mitochondrial inner membrane']

#
# membrane_size_specific = membrane_size_tr0[organelle1_m]
membrane_size_specific = volume_size_tr0[organelle1_m]
membrane_size_specific["mitochondrial outer membrane"] = membrane_size_specific["mitochondrial outer membrane"]/volume_size_tr0["mitochondrion"]
membrane_size_specific["mitochondrial inner membrane"] = membrane_size_specific["mitochondrial inner membrane"]/volume_size_tr0["mitochondrion"]

membrane_size_specific["sample_ID"] = list(membrane_size_specific.index)
membrane_size_specific["dilution rate (/h)"] = singleMapping(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["kinetic"], membrane_size_specific["sample_ID"])

x0 = "dilution rate (/h)"
for y0 in organelle1_m:
    title0 = 'result/figure/rose_membrane_pro_to_volume ' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=membrane_size_specific, lowess=True, height=4, aspect=1)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " pro volume in its organelle",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')















