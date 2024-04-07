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
pro_abundance = pro_abundance.rename(columns={'Systematic Name': 'gene'})
pro_abundance = pro_abundance.drop('Standard Name', axis=1)


pro_abundance1 = pd.read_excel("data/ProteomicsData_cell_systems_2018/mmc9_normal_plus_stress_conditions.xlsx")
pro_abundance1 = pro_abundance1.rename(columns={'Systematic Name': 'gene'})
pro_abundance1 = pro_abundance1.drop('Standard Name', axis=1)

######################################################
# input data from cell reports 2017
# note this data is obtained under exponential growth phases
pro_abundance2 = pd.read_excel("data/proteomics/protein_copy_cell_report_2017.xlsx")
# filter out one sample with very few total protein copy number
pro_abundance2 =pro_abundance2[[x for x in pro_abundance2.columns if x !="Chong et al. 2015 - Copy "]]


######################################################
# input data from zinc proteome.xlsx
pro_abundance3 = pd.read_excel("data/proteomics/zinc proteome.xlsx")
# filter out one sample with very few total protein copy number
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
# ID mapping between uniprot ID and gene locus IDs
[x for x in pro_abundance3['Uniprot ID'].tolist() if "-" in x]
pro_abundance3['Uniprot ID'] = pro_abundance3['Uniprot ID'].str.replace("-2","").str.replace("-3", "")
pro_abundance3['gene'] = multiMapping(id_mapping['GeneName'], id_mapping['Entry'], pro_abundance3['Uniprot ID'])
pro_abundance3 = pro_abundance3.drop('Uniprot ID', axis=1)



# a functiion to combine the different proteomics
def combineAbosluteAbundance(omics_combine_base, omics_new, remove_column="gene"):
    df_combine_auto = pd.merge(left=omics_combine_base, right=omics_new, left_on=['all_gene'], right_on=['gene'], how="left")
    # remove the duplicated
    df_combine_auto.pop(remove_column)
    return df_combine_auto
# combine data from different source?
# get all genes
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
# all_gene = set(omics_carl['gene'].tolist()) | set(omics_francesca['genes'].tolist()) | set(omics_tao1['gene'].tolist()) | set(omics_johan['gene'].tolist()) | set(omics_Tyler['gene'].tolist()) | set(omics_Rahul['gene'].tolist())
# all_gene = list(set(all_gene))
all_gene = sce_gene["GeneName"].tolist()
new_df = pd.DataFrame({"all_gene": all_gene})
new_df = new_df.dropna()
df_combine0 = combineAbosluteAbundance(new_df, pro_abundance, remove_column="gene")
df_combine1 = combineAbosluteAbundance(df_combine0, pro_abundance1, remove_column="gene")
df_combine2 = combineAbosluteAbundance(df_combine1, pro_abundance2, remove_column="gene")
df_combine3 = combineAbosluteAbundance(df_combine2, pro_abundance3, remove_column="gene")
df_combine3 = df_combine3.rename(columns={'all_gene': 'gene'})


# change the data as mass fraction for each protein per gram of total protein mass
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# change the unit from copy/cell or mmol/gDCW into g/g total protein?
df_combine3["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], df_combine3["gene"])
df_combine3 = df_combine3[~df_combine3["MW_Kda"].isna()]
all_colum = df_combine3.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']
omics_combine_input2 = df_combine3[all_colum2]
for x in all_colum2:
    omics_combine_input2[x] = df_combine3[x]*df_combine3["MW_Kda"]
omics_combine_input_mass_fraction = omics_combine_input2
for x in all_colum2:
    omics_combine_input_mass_fraction[x] = omics_combine_input2[x]/omics_combine_input2[x].sum()
omics_combine_input_mass_fraction["gene"] = df_combine3["gene"]
new_column = ["gene"] + all_colum2
mass_fraction_final = omics_combine_input_mass_fraction[new_column]
mass_fraction_final.to_excel("data/proteomics/mass_fraction_from_protein_copy.xlsx")



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


# # calculate the mass ratio
out = ProMassRatio_Organelle(mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_from_protein_copy.xlsx")


# calculate the volume ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_from_protein_copy.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_from_protein_copy.xlsx")
