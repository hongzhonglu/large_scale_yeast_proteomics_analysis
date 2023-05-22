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







# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]
volume_size_rosemery["sample_ID"] = list(volume_size_rosemery.index)
volume_size_rosemery["total_pro_volume"] = singleMapping(total_pro_volume['total_pro_volume'], total_pro_volume['sampleID'], volume_size_rosemery["sample_ID"])
column0 = list(volume_size_rosemery.columns)[0:139] + ["total_pro_volume"]
volume_size_rosemery_ratio = volume_size_rosemery[column0]
volume_size_rosemery_ratio1 = volume_size_rosemery_ratio.copy()

column1 = list(volume_size_rosemery.columns)[0:139]
for x in column1:
    print(x)
    volume_size_rosemery_ratio1[x] = volume_size_rosemery_ratio[x] / volume_size_rosemery_ratio["total_pro_volume"]





volume_size_rosemery_ratio1['sample_ID'] = list(volume_size_rosemery_ratio1.index)




# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery_ratio1, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation


# check other volume
# here we only explore the condition with only NH4 limitation
combine_data2 = combine_data[combine_data["Nitrogen source"] =="NH4"]
combine_data2 = combine_data2[combine_data2["limiting nutrient"] =="N"]
combine_data2.to_excel("data/proteomics/organell_protein_ratio_rosemary.xlsx")


# first calculate correlation coefficient
from scipy.stats import pearsonr
correlation = []
x0 = "dilution rate (/h)"
for y0 in column1:
    corr, _ = pearsonr(combine_data2[x0], combine_data2[y0])
    correlation.append(corr)
new_df = pd.DataFrame({"organelle":column1,"coefficient":correlation})
pd_null = new_df.sort_values(by=['coefficient'], ascending=False)
pd_null["coefficient_abs"] = abs(pd_null["coefficient"])
pd_null1 = pd_null[pd_null["coefficient_abs"] >= 0.85]

# bar plot
x0='organelle'
y0='coefficient'
plt.figure(figsize=(8, 6))
#sns.set_style('darkgrid')
sns.barplot(x=x0, y=y0, data=pd_null1, capsize=.2)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/organelle protein volume correlation with growth.pdf", bbox_inches='tight')






column_select = list(combine_data2.columns)
column_select1 = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet',
                  'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'cytosolic ribosome','mitochondrial ribosome','nucleolus']

# note: The nucleolus is a region found within the cell nucleus that is concerned with producing and assembling the cell's ribosomes.
# 'mitochondrial ribosome'


x0 = "dilution rate (/h)"
for y0 in column_select1:
    title0 = 'result/figure/rose_miu_ratio_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    #sns.lmplot(x=x0, y=y0, data=combine_data2,
    #           lowess=True,height=4, aspect=1)
    plt.figure(figsize=(4, 4))
    sns.lineplot(x=x0, y=y0, data=combine_data2, marker="o")
    plt.xlabel(x0, fontsize=15)
    plt.ylabel(y0, fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.savefig(title0,  bbox_inches='tight')

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







# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]
volume_size_rosemery["sample_ID"] = list(volume_size_rosemery.index)
volume_size_rosemery["total_pro_volume"] = singleMapping(total_pro_volume['total_pro_volume'], total_pro_volume['sampleID'], volume_size_rosemery["sample_ID"])
column0 = list(volume_size_rosemery.columns)[0:139] + ["total_pro_volume"]
volume_size_rosemery_ratio = volume_size_rosemery[column0]
volume_size_rosemery_ratio1 = volume_size_rosemery_ratio.copy()

column1 = list(volume_size_rosemery.columns)[0:139]
for x in column1:
    print(x)
    volume_size_rosemery_ratio1[x] = volume_size_rosemery_ratio[x] / volume_size_rosemery_ratio["total_pro_volume"]





volume_size_rosemery_ratio1['sample_ID'] = list(volume_size_rosemery_ratio1.index)




# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery_ratio1, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation


# check other volume
# here we only explore the condition with only NH4 limitation
combine_data2 = combine_data[combine_data["Nitrogen source"] =="NH4"]
combine_data2 = combine_data2[combine_data2["limiting nutrient"] =="N"]
combine_data2.to_excel("data/proteomics/ratio_of_organelle_volume_to_total_protein_volume.xlsx")

# first calculate correlation coefficient
from scipy.stats import pearsonr
correlation = []
x0 = "dilution rate (/h)"
for y0 in column1:
    corr, _ = pearsonr(combine_data2[x0], combine_data2[y0])
    correlation.append(corr)
new_df = pd.DataFrame({"organelle":column1,"coefficient":correlation})
pd_null = new_df.sort_values(by=['coefficient'], ascending=False)
pd_null["coefficient_abs"] = abs(pd_null["coefficient"])
pd_null1 = pd_null[pd_null["coefficient_abs"] >= 0.85]

# bar plot
x0='organelle'
y0='coefficient'
plt.figure(figsize=(8, 6))
sns.set_style('darkgrid')
sns.barplot(x=x0, y=y0, data=pd_null1, capsize=.2)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/organelle protein volume correlation with growth.pdf", bbox_inches='tight')






column_select = list(combine_data2.columns)
column_select1 = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet',
                  'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'cytosolic ribosome','mitochondrial ribosome','nucleolus']

# note: The nucleolus is a region found within the cell nucleus that is concerned with producing and assembling the cell's ribosomes.
# 'mitochondrial ribosome'


x0 = "dilution rate (/h)"
for y0 in column_select1:
    title0 = 'result/figure/rose_miu_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True,height=4, aspect=1)
    plt.xlabel(x0, fontsize=15)
    plt.ylabel(y0, fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.savefig(title0)


# add the yeast cell size data
# check the correlation between organelle proteins volume with cell size
df_curated = calculateCurationCoefficent()
df_curated["total_pro_volume"] = total_pro_volume["total_pro_volume"]
df_curated["total_pro_volume/cell_size"] = df_curated["total_pro_volume"]/df_curated["cell_size"]
# plot
x0="growth_rate"
y0="total_pro_volume/cell_size"
sns.lmplot(x=x0, y=y0, data=df_curated, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.axvline(x=0.18, color='k', linestyle='--')



"""
x0 = 'qGlucose (mmol/gDW h)'
for y0 in column_select1:
    title0 = 'result/figure/rose_qs_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True, height=4, aspect=1)
    plt.xlabel(x0, fontsize=15)
    plt.ylabel(y0, fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.axvline(x=3.2, color='k', linestyle='--')
    plt.savefig(title0)
"""




# add the yeast cell size data
# check the correlation between organelle proteins volume with cell size
df_curated = calculateCurationCoefficent()
df_curated["total_pro_volume"] = total_pro_volume["total_pro_volume"]
df_curated["total_pro_volume/cell_size"] = df_curated["total_pro_volume"]/df_curated["cell_size"]
# plot
x0="growth_rate"
y0="total_pro_volume/cell_size"
sns.lmplot(x=x0, y=y0, data=df_curated, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.axvline(x=0.18, color='k', linestyle='--')


