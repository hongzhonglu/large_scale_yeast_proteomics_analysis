# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns


# input the membrane size data
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_tao2.xlsx") # curated based on cell size under different growth rate


volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])
volume_size_tr0 = volume_size_tr0.iloc[2: ,:]



# calculate the total volume of proteins
protein_copy_all1 = pd.read_excel("data/proteomics/protein_copy_tao_nc.xlsx")
# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

protein_copy_all1["pro_volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],protein_copy_all1['gene'])


sample_ID = list(protein_copy_all1.columns)
sample_ID = [x for x in sample_ID if x != "gene"]
sample_ID = [x for x in sample_ID if x != "pro_volume"]

volume_list = []
for x in sample_ID:
    ss1 = protein_copy_all1[[x,"pro_volume"]]
    ss1["value"] = ss1[x]*ss1["pro_volume"]
    ss1 = ss1[~ss1["value"].isna()]
    sum0 = sum(ss1['value'])
    # change nm^3 into um^3
    total_volume_um = sum0 / 1e9
    volume_list.append(total_volume_um)

# creat a new dataframe
total_pro_volume = pd.DataFrame({"sampleID":sample_ID,"total_pro_volume":volume_list})










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
combine_data2["sample_ID"] = [5,5,115,115,30,30,50,50]


# add the total protein volume
combine_data2["total_pro_volume"] = total_pro_volume["total_pro_volume"].tolist()
combine_data2["sample"] = total_pro_volume["sampleID"].tolist()

volume_size_ratio1 = combine_data2.copy()

column1 = list(combine_data2.columns[0:13])
for x in column1:
    print(x)
    volume_size_ratio1[x] = combine_data2[x] / combine_data2["total_pro_volume"]


# plot the figure in ratio
x0 = "sample_ID"
for y0 in column_select1:
    title0 = 'result/figure/tao2_' + y0 + '.pdf'
    print(title0)
    plt.figure()
    sns.barplot(x=x0, y=y0, data=combine_data2, capsize=.2)
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " protein volume",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(rotation=90)
    plt.savefig(title0, bbox_inches='tight')






# plot the figure in ratio
x0 = "sample_ID"
for y0 in column_select1:
    title0 = 'result/figure/tao2_ratio_' + y0 + '.pdf'
    print(title0)
    plt.figure()
    sns.barplot(x=x0, y=y0, data=volume_size_ratio1, capsize=.2)
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " per total protein volume",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xticks(rotation=90)
    plt.savefig(title0, bbox_inches='tight')









