# Note:
# which protein on the membrane
# 1) all the transporter protein
# 2) all the membrane-related protein


# import self function
from src.protein_process import *

# input the protein abundance data
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx")
# test
omics_combine_auto = omics_combine_auto[['all_gene'] + [x for x in omics_combine_auto.columns if "_M" in x]]
omics_combine_auto.columns = omics_combine_auto.columns.str.replace('all_gene', 'gene')
protein_copy1 = omics_combine_auto.copy()





# change the data as mass fraction for each protein per gram of total protein mass
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# change the unit from copy/cell or mmol/gDCW into g/g total protein?
protein_copy1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], protein_copy1["gene"])
protein_copy1 = protein_copy1[~protein_copy1["MW_Kda"].isna()]
all_colum = protein_copy1.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']
omics_combine_input2 = protein_copy1[all_colum2]
for x in all_colum2:
    omics_combine_input2[x] = protein_copy1[x]*protein_copy1["MW_Kda"]
omics_combine_input_mass_fraction = omics_combine_input2
for x in all_colum2:
    omics_combine_input_mass_fraction[x] = omics_combine_input2[x]/omics_combine_input2[x].sum()
omics_combine_input_mass_fraction["gene"] = protein_copy1["gene"]
new_column = ["gene"] + all_colum2
mass_fraction_final = omics_combine_input_mass_fraction[new_column]


# how to further calculation the protein volume ratio and protein area ratio of main organelle
# change the unit from g/g into mol/g?
mass_fraction = mass_fraction_final.copy()
mass_fraction["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], mass_fraction["gene"])
all_colum = mass_fraction.columns
all_colum2 = [x for x in all_colum1 if x !='gene']
protein_in_mol = mass_fraction[all_colum2]
for x in all_colum2:
    protein_in_mol[x] = 1000 * protein_in_mol[x] / mass_fraction["MW_Kda"]
protein_in_mol["gene"] = mass_fraction["gene"]
new_column = ["gene"] + all_colum2
protein_in_mol = protein_in_mol[new_column]

# # calculate the mass ratio
out = ProMassRatio_Organelle(mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_jianye_test2.xlsx")


# calculate the volume ratio
s2 =Pro_3D_Volume_Ratio_Cal(protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_jianye_test2.xlsx")


# calculate the membrane ratio
s2 = Pro_Membrance_Ratio_Cal(protein_in_mol)
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_jianye_test2.xlsx")

