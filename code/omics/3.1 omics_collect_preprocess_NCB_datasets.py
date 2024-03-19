# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

import sys

# import self function
from src.mainFunction import *
from src.protein_process import *


# some general datasets
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

# absolute part
proteomics_NCB1 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_sce.xlsx", sheet_name="Table 10a. abs_prot_SC_CENPK")
proteomics_NCB2 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_sce.xlsx", sheet_name="Table 10b. abs_prot_SC_FY4")

proteomics_NCB1 = proteomics_NCB1[["Entry","geneID","mean"]]
proteomics_NCB1.columns = ["Entry","geneID","sce_CEN.PK_batch"]
proteomics_NCB2 = proteomics_NCB2[["Entry","geneID","mean"]]
proteomics_NCB2.columns = ["Entry","geneID","sce_FY4_batch_miu=0.255"]

# 1 relative data, change it as the absolute value
proteomics_NCB3 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_sce.xlsx", sheet_name="Table 11a. rel_prot_SC_FY4") # sce_FY4
proteomics_NCB4 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_sce.xlsx", sheet_name="Table 11c. rel_prot_SC_resp_inh") # sce_cenpk
proteomics_NCB3 = proteomics_NCB3.iloc[:,0:14] # sce_FY4
proteomics_NCB4 = proteomics_NCB4.iloc[:,0:9] # sce_cenpk

proteomics_NCB3["ref_abs"] = multiMapping(proteomics_NCB2["sce_FY4_batch_miu=0.255"],proteomics_NCB2["Entry"],proteomics_NCB3["Entry"])
proteomics_NCB3 = proteomics_NCB3.replace(to_replace='None', value=np.nan).dropna()
proteomics_NCB3['ref_abs'] = proteomics_NCB3['ref_abs'].apply(lambda x: float(x))
target_column_FY4 = list(proteomics_NCB3.columns)[2:14]
proteomics_NCB30 = proteomics_NCB3[["Entry"]+target_column_FY4]
for xx in target_column_FY4:
    print(xx)
    proteomics_NCB30[xx] = proteomics_NCB3.ref_abs * proteomics_NCB3[xx]
# unify the name
target_column_FY40 = ["sce_FY4_" + x for x in target_column_FY4]
miu_measured = ["miu=0.08", "miu=0.16", "miu=0.22", "miu=0.28"]*3
target_column_FY40 = [x+"_"+y for x, y in zip(target_column_FY40, miu_measured)]


proteomics_NCB30.columns = ["Entry"] + target_column_FY40
# merge all the dataset
df_combine1 = pd.merge(left=proteomics_NCB2, right=proteomics_NCB30, left_on=['Entry'], right_on=['Entry'], how="left")



# 2 relative data, change it as the absolute value
proteomics_NCB4 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_sce.xlsx", sheet_name="Table 11c. rel_prot_SC_resp_inh") # sce_cenpk
proteomics_NCB4 = proteomics_NCB4.iloc[:,0:9] # sce_cenpk

proteomics_NCB4["ref_abs"] = multiMapping(proteomics_NCB1["sce_CEN.PK_batch"],proteomics_NCB1["Entry"],proteomics_NCB4["Entry"])
proteomics_NCB4 = proteomics_NCB4.replace(to_replace='None', value=np.nan).dropna()
proteomics_NCB4['ref_abs'] = proteomics_NCB4['ref_abs'].apply(lambda x: float(x))
target_column_CEN = list(proteomics_NCB4.columns)[1:9]
proteomics_NCB40 = proteomics_NCB4[["Entry"]+target_column_CEN]
for xx in target_column_CEN:
    print(xx)
    proteomics_NCB40[xx] = proteomics_NCB4.ref_abs * proteomics_NCB4[xx]
# unify the name
target_column_CEN0 = ["sce_CEN.PK_" + x for x in target_column_CEN]
proteomics_NCB40.columns = ["Entry"] + target_column_CEN0
# merge the dataset
df_combine2 = pd.merge(left=proteomics_NCB1, right=proteomics_NCB40, left_on=['Entry'], right_on=['Entry'], how="left")
df_combine2 = df_combine2.rename(columns={'sce_CEN.PK_batch': 'sce_CEN.PK_batch_miu=0.391'})

df_combine1.to_excel("data/mass_fraction_NCB.xlsx")
df_combine2.to_excel("data/mass_fraction_NCB2.xlsx")

# merge all the dataset
mass_fraction_NCB = pd.merge(left=df_combine1, right=df_combine2, left_on=['Entry'], right_on=['Entry'], how="outer")


mass_fraction_NCB = mass_fraction_NCB.drop(['Entry', 'geneID_y'], axis=1)
mass_fraction_NCB = mass_fraction_NCB.rename(columns={'geneID_x': 'gene'})
mass_fraction_final = mass_fraction_NCB.copy()
mass_fraction_final.to_excel("data/proteomics/mass_fraction_NCB.xlsx")

# the above unit is g specific protein/ g total protein. However, the mass percentage of protein in some condition was not measured.
# change the unit as protein abundance data in the unit of mmol/g DCW ??
# calculate mass ratio of protein from each organelle per total protein mass



# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000

# how to further calculation the protein volume ratio and protein area ratio of main organelle
# change the unit from g/g into mol/g?
mass_fraction = mass_fraction_final.copy()
mass_fraction["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], mass_fraction["gene"])
all_colum = mass_fraction.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']
protein_in_mol = mass_fraction[all_colum2]
for x in all_colum2:
    protein_in_mol[x] = 1000 * protein_in_mol[x] / mass_fraction["MW_Kda"]
protein_in_mol["gene"] = mass_fraction["gene"]
new_column = ["gene"] + all_colum2
protein_in_mol = protein_in_mol[new_column]


# test the above code
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_NCB.xlsx")

# calculate the volume ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_NCB.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_NCB.xlsx")



