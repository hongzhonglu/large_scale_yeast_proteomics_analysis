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
batch_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/batch_ExperimentalData.csv")
batch1 = batch_ExperimentalData.iloc[0,].tolist()
batch2 = list(batch_ExperimentalData.columns)
batch3 = [str(x) + "_miu=" + str(y)for x, y in zip(batch2, batch1)]
batch3[0] = "gene"
batch_ExperimentalData.columns = batch3
batch_ExperimentalData0 = batch_ExperimentalData.iloc[1:,]


translation_inhibit_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/TI_ExperimentalData.csv")
translation_inhibit1 = translation_inhibit_ExperimentalData.iloc[0,].tolist()
translation_inhibit2 = list(translation_inhibit_ExperimentalData.columns)
translation_inhibit3 = [str(x) + "_miu=" + str(y)for x, y in zip(translation_inhibit2, translation_inhibit1)]
translation_inhibit3[0] = "gene"
translation_inhibit_ExperimentalData.columns = translation_inhibit3
translation_inhibit_ExperimentalData0 = translation_inhibit_ExperimentalData.iloc[1:,]


chemostat_inhibit_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/chemostat_ExperimentalData.csv")
chemostat_inhibit1 = chemostat_inhibit_ExperimentalData.iloc[0,].tolist()
chemostat_inhibit2 = list(chemostat_inhibit_ExperimentalData.columns)
chemostat_inhibit3 = [str(x) + "_miu=" + str(y)for x, y in zip(chemostat_inhibit2, chemostat_inhibit1)]
chemostat_inhibit3[0] = "gene"
chemostat_inhibit_ExperimentalData.columns = chemostat_inhibit3
chemostat_inhibit_ExperimentalData0 = chemostat_inhibit_ExperimentalData.iloc[1:,]



# combine all the above dataset together
mass_fraction_ibrahim = pd.merge(left=batch_ExperimentalData0, right=translation_inhibit_ExperimentalData0, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_ibrahim = pd.merge(left=mass_fraction_ibrahim, right=chemostat_inhibit_ExperimentalData0, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_final = mass_fraction_ibrahim.copy()




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






#
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_ibrahim.xlsx")

# calculate the volume ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_jibrahim.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_ibrahim.xlsx")




