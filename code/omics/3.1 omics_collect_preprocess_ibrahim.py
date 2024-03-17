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

# test the above code
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_ibrahim, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_ibrahim.xlsx")
