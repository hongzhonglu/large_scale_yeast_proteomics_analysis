# Note:
# which protein on the membrane
# 1) all the transporter protein
# 2) all the membrane-related protein


# import self function
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")


# only analyze the rosemary datasets under NH4 limitation?
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all_rosemary = protein_copy_all1[Sample_ID_select+["gene"]]
# curate the protein copies based on newly calculated size
df_curated = calculateCurationCoefficent()
curation_coefficent = df_curated["curation_coefficent"].to_list()
for i, sid in enumerate(Sample_ID_select):
    print(i, sid, curation_coefficent[i])
    protein_copy_all_rosemary[sid] = protein_copy_all_rosemary[sid]*curation_coefficent[i]

protein_copy_all_rosemary.to_excel("data/proteomics/all_protein_copy_rosemary.xlsx")

s2 = ProMembraneCal(protein_copy_all_rosemary)
s2.to_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v3.xlsx")


# jianye datasets
all_columns = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in all_columns if "_M" in x]
protein_copy_jianye = protein_copy_all1[Sample_ID_select+["gene"]]
s2 = ProMembraneCal(protein_copy_jianye)
s2.to_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")


# all datasets
s2 =ProMembraneCal(protein_copy_all1)
s2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")


