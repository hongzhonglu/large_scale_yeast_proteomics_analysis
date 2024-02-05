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
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])




# calculate the total volume of proteins
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy_rosemary.xlsx")
# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

protein_copy_all1["pro_volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],protein_copy_all1['gene'])
sample_ID = list(protein_copy_all1.columns)
sample_ID = sample_ID[1:19]
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



nucleolus_pro_volume = protein_copy_all1.copy()
for x in sample_ID:
    ss1 = protein_copy_all1[[x,"pro_volume"]]
    ss1["value"] = ss1[x]*ss1["pro_volume"] / 1e9
    nucleolus_pro_volume[x] = ss1["value"]


volume_size_tr = nucleolus_pro_volume.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[19])
volume_size_tr0 = volume_size_tr0.iloc[1:19]

# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]
volume_size_rosemery["sample_ID"] = list(volume_size_rosemery.index)
volume_size_rosemery["total_pro_volume"] = singleMapping(total_pro_volume['total_pro_volume'], total_pro_volume['sampleID'], volume_size_rosemery["sample_ID"])
column0 = list(volume_size_rosemery.columns)[0:6860] + ["total_pro_volume"]
volume_size_rosemery_ratio = volume_size_rosemery[column0]
volume_size_rosemery_ratio1 = volume_size_rosemery_ratio.copy()

volume_size_rosemery = volume_size_rosemery.dropna(axis=1)
column1 = list(volume_size_rosemery.columns)[0:3128]
for x in column1:
    print(x)
    volume_size_rosemery_ratio1[x] = volume_size_rosemery_ratio[x] / volume_size_rosemery_ratio["total_pro_volume"]





volume_size_rosemery_ratio1['sample_ID'] = list(volume_size_rosemery_ratio1.index)

# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery_ratio1, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")

# further analysis
# first calculate correlation coefficient
from scipy.stats import pearsonr
correlation = []
x0 = "dilution rate (/h)"
for y0 in column1:
    corr, _ = pearsonr(combine_data[x0], combine_data[y0])
    correlation.append(corr)
new_df = pd.DataFrame({"gene":column1,"coefficient":correlation})
pd_null = new_df.sort_values(by=['coefficient'], ascending=False)

# focus on the genes with high coefficient
pd_null_filter1 = pd_null[pd_null['coefficient'] >= 0.9]
pd_null_filter2 = pd_null[pd_null['coefficient'] <= -0.8]

gene_list1 = ','.join(pd_null_filter1['gene'].tolist())
print(gene_list1)
gene_list2 = ','.join(pd_null_filter2['gene'].tolist())
print(gene_list2)