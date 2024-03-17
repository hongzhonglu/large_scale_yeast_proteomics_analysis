import matplotlib.pyplot as plt
import os


# import self function
from src.protein_process import *

# Part 2 Collect all the data in the molecular/cell
# Generally, there are three sources, SGD, cell system and another paper. The SGD data use the median value from cell system.
# input data from SGD
# pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')

# input data from cell system, 2018
pro_abundance = pd.read_excel("data/ProteomicsData_cell_systems_2018/mmc5_normal_conditions.xlsx")
col0 = list(pro_abundance.columns)
col0 = [x for x in col0 if x != "Coefficient of Variation"]
pro_abundance = pro_abundance[col0]
col0 = col0[0:4] + [x+"_ref" for x in col0[4:]]
pro_abundance.columns = col0

pro_abundance1 = pd.read_excel("data/ProteomicsData_cell_systems_2018/mmc9_normal_plus_stress_conditions.xlsx")


# input data from cell reports 2017
# note this data is obtained under exponential growth phases
pro_abundance2 = pd.read_excel("data/proteomics/protein_copy_cell_report_2017.xlsx")
# filter out one sample with very few total protein copy number
pro_abundance2 =pro_abundance2[[x for x in pro_abundance2.columns if x !="Chong et al. 2015 - Copy "]]



# combine data from different source?
protein_copy = pd.merge(left=pro_abundance, right=pro_abundance1, left_on=['Systematic Name'], right_on=['Systematic Name'], how="outer")
protein_copy1 = pd.merge(left=protein_copy, right=pro_abundance2, left_on=['Systematic Name'], right_on=['gene'], how="outer")
col2 = list(protein_copy1.columns)
col2 = [x for x in col2 if "_x" not in x]
col2 = [x for x in col2 if "_y" not in x]
col2 = [x for x in col2 if x !="gene"]
protein_copy1 = protein_copy1[col2]
protein_copy1 = protein_copy1.rename(columns={'Systematic Name': 'gene'})

protein_copy1.to_excel("data/proteomics/protein_copy_combine.xlsx",index=False)






# change the data as mass fraction for each protein per gram of total protein mass
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")


protein_copy1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], protein_copy1["gene"])

protein_copy1 = protein_copy1[~protein_copy1["MW_Kda"].isna()]

all_colum = protein_copy1.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']

omics_combine_input2 = protein_copy1[all_colum2]


# change copy/cell into g/g total protein:
for x in all_colum2:
    omics_combine_input2[x] = protein_copy1[x]*protein_copy1["MW_Kda"]
omics_combine_input_mass_fraction = omics_combine_input2

for x in all_colum2:
    omics_combine_input_mass_fraction[x] = omics_combine_input2[x]/omics_combine_input2[x].sum()

omics_combine_input_mass_fraction["gene"] = protein_copy1["gene"]
new_column = ["gene"] + all_colum2
mass_fraction_cell_system_2018 = omics_combine_input_mass_fraction[new_column]


# how to further calculation the protein volume ratio and protein area ratio of main organelle
# change the unit from g/g into mol/g?
mass_fraction_cell_system_2018["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], mass_fraction_cell_system_2018["gene"])

all_colum = mass_fraction_cell_system_2018.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']
protein_in_mol = mass_fraction_cell_system_2018[all_colum2]
# change copy/cell into g/g total protein:
for x in all_colum2:
    protein_in_mol[x] = 1000 * protein_in_mol[x] / mass_fraction_cell_system_2018["MW_Kda"]

protein_in_mol["gene"] = mass_fraction_cell_system_2018["gene"]
new_column = ["gene"] + all_colum2

protein_in_mol = protein_in_mol[new_column]

# test the above code
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_cell_system_2018, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_cell_system_2018.xlsx")



# calculate the membrane ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_copy=protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_cell_system_2018.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_cell_system_2018.xlsx")

