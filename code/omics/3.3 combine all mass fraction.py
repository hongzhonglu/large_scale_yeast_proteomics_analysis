import matplotlib.pyplot as plt
import os


# import self function
from src.protein_process import *

# input the new absolute proteomics
mass_fraction_NCB = pd.read_excel("data/proteomics/mass_fraction_NCB.xlsx")
mass_fraction_NCB = mass_fraction_NCB.iloc[:,1:]

mass_fraction_ibrahim = pd.read_excel("data/proteomics/mass_fraction_ibrahim.xlsx")
mass_fraction_ibrahim = mass_fraction_ibrahim.iloc[:,1:]

mass_fraction_from_protein_copy = pd.read_excel("data/proteomics/mass_fraction_from_protein_copy.xlsx")
mass_fraction_from_protein_copy = mass_fraction_from_protein_copy.iloc[:,1:]

mass_fraction_others = pd.read_excel("data/proteomics/mass_fraction_others.xlsx")
mass_fraction_others = mass_fraction_others.iloc[:,1:]

mass_fraction_all = pd.merge(left=mass_fraction_others, right=mass_fraction_ibrahim, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_all = pd.merge(left=mass_fraction_all, right=mass_fraction_NCB, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_all = pd.merge(left=mass_fraction_all, right=mass_fraction_from_protein_copy, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_final = mass_fraction_all.copy()

# Save
mass_fraction_final.to_excel("data/proteomics/mass_fraction_combine.xlsx", index=False) # the unit the mmol/gDW



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


# calculate the mass ratio
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_combine_test.xlsx")


# calculate the volume ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_in_mol, compartment_type="organelle")
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_combine.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_combine.xlsx")


