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
# input the membrane size data
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment.xlsx")
volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])




# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]



volume_size_rosemery["sample_ID"] = list(volume_size_rosemery.index)
# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation


# check other volume
# here we only explore the condition with only NH4 limitation
combine_data2 = combine_data[combine_data["Nitrogen source"] =="NH4"]
combine_data2 = combine_data2[combine_data2["limiting nutrient"] =="N"]

column_select = list(combine_data2.columns)
column_select1 = [x for x in column_select if "membrane" not in x]
column_select1 = column_select1[0:105]

x0 = "dilution rate (/h)"
for y0 in column_select1:
    title0 = 'result/figure/rose_miu_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True)
    plt.xlabel(x0)
    plt.ylabel(y0)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.savefig(title0)





x0 = 'qO2 (mmol/gDW h)'
for y0 in column_select1:
    title0 = 'result/figure/rose_qo2_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True)
    plt.xlabel(x0)
    plt.ylabel(y0)
    plt.axvline(x=5.4, color='k', linestyle='--')
    plt.savefig(title0)



x0 = 'qGlucose (mmol/gDW h)'
for y0 in column_select1:
    title0 = 'result/figure/rose_qs_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True)
    plt.xlabel(x0)
    plt.ylabel(y0)
    plt.axvline(x=3.2, color='k', linestyle='--')
    plt.savefig(title0)






# only take Jianye physiology dataset
physiology_Jianye = physiology_data[physiology_data["kinetic"].str.contains("_M")]
# only take Jianye proteomics
volume_size_Jianye = volume_size_tr0[volume_size_tr0.index.str.contains("_M")]

volume_size_Jianye["sample_ID"] = list(volume_size_Jianye.index)
# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_Jianye, right=physiology_Jianye, left_on=['sample_ID'], right_on=['kinetic'], how="left")
combine_data2 = combine_data

column_select = list(combine_data2.columns)
column_select1 = [x for x in column_select if "membrane" not in x]
column_select1 = column_select1[0:105]


x0 = "dilution rate (/h)"
for y0 in column_select1:
    title0 = 'result/figure/jianye_miu_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2,
               lowess=True)
    plt.xlabel(x0)
    plt.ylabel(y0)
    plt.axvline(x=0.3, color='k', linestyle='--')
    plt.savefig(title0)