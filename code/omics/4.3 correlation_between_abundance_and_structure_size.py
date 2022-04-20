# check the following correlation between:
# 1. organelle volume, phenotype data
# 2. organelle membrane, phenotype data
# 3. organelle membrane, the related protein abundance


import matplotlib.pyplot as plt
import os
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *
import seaborn as sns
from scipy.stats import pearsonr


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



# check the correlation between protein abundance and sectional surface size
abundance_organelle = pd.read_excel("data/proteomics/protein_abundance_across_compartment.xlsx")

Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
abundance_organelle = abundance_organelle[Sample_ID_select +["compartment"]]
abundance_organelle = abundance_organelle[abundance_organelle["compartment"].isin(column_select)]
membrane_size0 = membrane_size[Sample_ID_select +["compartment"]]
membrane_size0 = membrane_size0[membrane_size0["compartment"].isin(column_select)]

correlation = []
for xx in column_select:
    print(xx)
    ss1 = membrane_size0[membrane_size0["compartment"]==xx].values.tolist()[0][0:18]
    ss2 = abundance_organelle[abundance_organelle["compartment"]==xx].values.tolist()[0][0:18]
    pccs = pearsonr(ss1, ss2)
    correlation.append(pccs[0])

    title0 = 'result/figure/correlaton_' + xx + '.pdf'
    print(title0)
    #plt.figure()
    plt.figure(figsize=(5, 5))
    ss2 = [x/100000 for x in ss2]
    plt.scatter(ss2, ss1)
    plt.xlabel("total protein copy/10^5",fontsize=18)
    plt.ylabel("total protein sectional area",fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(xx, y=1.01)
    plt.savefig(title0, bbox_inches='tight')


