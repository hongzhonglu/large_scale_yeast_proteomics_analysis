# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
# input the membrane size data
# volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation.xlsx") # not curated
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_jianye_C_limitation.xlsx") # curated based on cell size under different growth rate


volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])
volume_size_tr0 = volume_size_tr0.iloc[2: ,:]

# calculate the total volume of proteins
protein_copy_all1 = pd.read_excel("data/proteomics/protein_copy_jianye.xlsx")

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

protein_copy_all1["pro_volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],protein_copy_all1['gene'])
sample_ID = list(protein_copy_all1.columns)

sample_ID = sample_ID[0:9]
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

volume_size_tr0["total_pro_volume"] = total_pro_volume["total_pro_volume"].tolist()



# only take jianye physiology dataset
physiology_jianye = physiology_data[physiology_data["kinetic"].str.contains("_M")]
# only take jianye proteomics
volume_size_jianye_ratio1 = volume_size_tr0.copy()

column1 = list(volume_size_tr0.columns)[0:140]

for x in column1:
    if x != "total_pro_volume":
        print(x)
        volume_size_jianye_ratio1[x] = volume_size_tr0[x] / volume_size_tr0["total_pro_volume"]

volume_size_jianye_ratio1['sample_ID'] = list(volume_size_jianye_ratio1.index)
volume_size_jianye_ratio1.to_excel("data/proteomics/organell_protein_ratio_jianye.xlsx")



# combine the physiological datasets and proteomics datasets
physiology_jianye = physiology_data[physiology_data["kinetic"].str.contains("_M")]
combine_data = pd.merge(left=volume_size_jianye_ratio1, right=physiology_jianye, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation


# check other volume
column_select = list(combine_data.columns)
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

combine_data = combine_data[["dilution rate (/h)"] + column_select1]
x0 = "dilution rate (/h)"
for y0 in column_select1:
    title0 = 'result/figure/jianye_miu_' + y0 + ' ratio per total protein volume.pdf'
    print(title0)
    #plt.figure()
    if y0 != "mitochondrial ribosome":
        sns.lmplot(x=x0, y=y0, data=combine_data,
                   lowess=True, height=4, aspect=1)
        plt.xlabel(x0, fontsize=15)
        plt.ylabel(y0, fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlim(0, 0.4)
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.savefig(title0,bbox_inches='tight')
    else:
        plt.figure(figsize=(4, 4))
        sns.lineplot(x=x0, y=y0, data=combine_data, marker="o")
        plt.xlabel(x0, fontsize=15)
        plt.ylabel(y0, fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlim(0, 0.4)
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.savefig(title0,bbox_inches='tight')



combine_data.to_excel("data/proteomics/organell protein volume ratio relative to total protein volume from Jianye.xlsx")

# try to plot all the result together
combine_data["sample_ID"] = volume_size_jianye_ratio1["sample_ID"].tolist()

combine_data2 = combine_data.transpose()
combine_data2.columns = volume_size_jianye_ratio1["sample_ID"].tolist()
combine_data2 = combine_data2.iloc[1:14,:]




# plot
df_in = combine_data2
x0='D=0.284_M'
y0='D=0.379_M'

df_in["Relative change"] = 2*(df_in[y0]-df_in[x0])/(df_in[x0]+df_in[y0])*100
df_in = df_in.sort_values(by=['Relative change'], ascending=False)
df_in["compartment"] = list(df_in.index)

plt.figure()
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=df_in)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.plot([0, 0.45], [0, 0.45], linewidth=2)
plt.xlim(0,0.45)
plt.ylim(0,0.45)
# annotate the dataset points
xs=df_in[x0].tolist()
ys=df_in[y0].tolist()
tlab=list(df_in.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()
plt.savefig("result/figure/organelle volume correlation analysis for " + y0 + " vs " + x0 + ".pdf", bbox_inches='tight')


# bar plot
x00='compartment'
y00='Relative change'
plt.figure()
sns.set_style('darkgrid')
sns.barplot(x=x00, y=y00, data=df_in, capsize=.2)
plt.xlabel(x00, fontsize=12)
plt.ylabel(y00 + "_" + y0 + " vs " + x0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/organelle volume relative change for " + y0 + " vs " + x0 + ".pdf", bbox_inches='tight')


