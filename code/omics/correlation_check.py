# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
# input the volume size data
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment.xlsx")
volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])

# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]


# only take NH4 limitation
physiology_rosemery2 = physiology_rosemery[physiology_rosemery["Nitrogen source"] =="NH4"]
physiology_rosemery2 = physiology_rosemery2[physiology_rosemery2["limiting nutrient"] =="N"]
column0 = list(physiology_rosemery2.columns)
column1 = [x for i,x in enumerate(column0) if i in [2,3,4,8,9,14,15,22,23]]
volume_size_rosemery2 = volume_size_rosemery[volume_size_rosemery.index.isin(physiology_rosemery2["kinetic"])]

size_column = list(volume_size_rosemery2.columns)
zeors_array = np.zeros( (len(column1), len(size_column)) )


from scipy.stats import pearsonr
for i in range(0,len(column1)):
    print(i)
    for j in range(0,len(size_column)):
        print(i, j)
        x = column1[i]
        y = size_column[j]
        pccs = pearsonr(physiology_rosemery2[x], volume_size_rosemery2[y])
        pearson_c = pccs[0]
        zeors_array[i,j] = pearson_c


df = pd.DataFrame(zeors_array, columns=size_column, index=column1)
df.to_excel("result/correlation_analysis_organelle_volume.xlsx")





# input the membrane size data
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment.xlsx")
membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])

# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
membrane_size_rosemery = membrane_size_tr0[membrane_size_tr0.index.str.contains("prot.")]


# only take NH4 limitation
physiology_rosemery2 = physiology_rosemery[physiology_rosemery["Nitrogen source"] =="NH4"]
physiology_rosemery2 = physiology_rosemery2[physiology_rosemery2["limiting nutrient"] =="N"]
column0 = list(physiology_rosemery2.columns)
column1 = [x for i,x in enumerate(column0) if i in [2,3,4,8,9,14,15,22,23]]
membrane_size_rosemery2 = membrane_size_rosemery[membrane_size_rosemery.index.isin(physiology_rosemery2["kinetic"])]

size_column = list(membrane_size_rosemery2.columns)
zeors_array = np.zeros( (len(column1), len(size_column)) )


from scipy.stats import pearsonr
for i in range(0,len(column1)):
    print(i)
    for j in range(0,len(size_column)):
        print(i, j)
        x = column1[i]
        y = size_column[j]
        pccs = pearsonr(physiology_rosemery2[x], membrane_size_rosemery2[y])
        pearson_c = pccs[0]
        zeors_array[i,j] = pearson_c


df = pd.DataFrame(zeors_array, columns=size_column, index=column1)

column_select =['mitochondrial outer membrane',
 'prospore membrane',
 'endoplasmic reticulum membrane',
 'mitochondrial inner membrane',
 'Golgi membrane',
 'mitochondrial membrane',
 'late endosome membrane',
 'peroxisomal membrane',
 'vacuolar membrane',
 'nuclear membrane',
 'endosome membrane',
 'nuclear inner membrane']

df = df[column_select]
df = df[~df.index.str.contains("Nitrogen")]
df0 = df.transpose()
ax = sns.heatmap(df0, cmap="YlGnBu")
plt.savefig("result/figure/correlation_analysis_membrane.pdf", bbox_inches='tight')



