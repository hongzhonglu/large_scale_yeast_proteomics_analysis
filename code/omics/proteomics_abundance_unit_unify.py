# this script is to rescale proteomics datasets
# 2021-11-16


import sys
sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
# import self function
from src.mainFunction import *



# coeffcient change

# then how to transform the absolute abundance into one cell?
# according to this paper https://www.pnas.org/content/107/3/999, we can get the density of yeast cell 1.1029 ± 0.0026 g/mL
# assume
Vgdw = 1.7 # unit is mL/gDW Ibrahim
# average volume of a cell
Vcell = 82 #µm3
Ncell = Vgdw/Vcell*1e12 #2.07e7
# set a coefficent to transform the mmol/g.Biomass into molecular/cell
# this step should be careful in the late calculation
coefficient0 = 6.022e20/Ncell






# part 1
# Input the datasets from paxDB
# update the unit
abundance = pd.read_csv("data/info_yeast/abundance_table.csv")
# Get the molecular weight data
mw = pd.read_csv("data/info_yeast/all_peptides_mw.csv")
mw["MW_Kda"] = mw["MW"]/1000
abundance["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance["genes"])
# scale the dataset
coefficient_protein_in_biomass = 0.45
abundance["abundance_w1"] = abundance["MW_Kda"]*abundance["abundance"]
w_total = sum(abundance["abundance_w1"])
abundance["abundance_scale"] = abundance["abundance_w1"]/w_total
abundance["abundance_per_biomass"] = abundance["abundance_scale"]*coefficient_protein_in_biomass # g/g biomass
abundance["abundance_per_biomass2"] = abundance["abundance_per_biomass"]/abundance["MW_Kda"]# #mmol/g biomass
abundance["copy_per_cell"] = abundance["abundance_per_biomass2"]*coefficient0 # Molecular/cell






# part2
# process the experimentally proteomic datasets
# This original dataset is from www.pnas.org/cgi/doi/10.1073/pnas.1918216117
abundance_ex = pd.read_excel("data/yeast_proteomics_example.xlsx")
abundance_ex = abundance_ex[["GeneNameOrdered","average(g/gDW)"]]
abundance_ex.columns = ["genes","absolute_abundance"]

# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000

abundance_ex["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex["genes"])

# Some species process. It could find that some rows with two proteins. Here assume the two proteins are identical in abundance，taking half of total protein abundance.
abudance_ex_1 = abundance_ex[~abundance_ex["MW_Kda"].isna()]

abudance_ex_check = abundance_ex[abundance_ex["MW_Kda"].isna()]
gene0 = []
abundance0 =[]
for i, x in abudance_ex_check.iterrows():
    print(i,x)
    s = x["genes"].split("; ")
    gene0 = gene0 + s
    v = [x["absolute_abundance"]/2]*2
    abundance0 = abundance0 + v
abudance_ex_2 = pd.DataFrame({"genes":gene0,"absolute_abundance":abundance0})
abudance_ex_2["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abudance_ex_2["genes"])

# combine the above two datasets
abundance_ex_corrected = pd.concat([abudance_ex_1, abudance_ex_2], axis=0)


# change the unit from g/g into mmol/g
abundance_ex_corrected["abundance_mmol"] = abundance_ex_corrected["absolute_abundance"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
one_strange = abundance_ex_corrected[abundance_ex_corrected['genes']=='YMR142C']
# here we should remove this protein in the following analysis.
# else, we should have an automatic way to filter out the abnormal values as we did not know which value is still not good!
abundance_ex_corrected = abundance_ex_corrected[~(abundance_ex_corrected['genes']=='YMR142C')]
abundance_ex_corrected["copy_per_cell"] = abundance_ex_corrected["abundance_mmol"] * coefficient0 # Molecular/cell
abundance_ex_corrected.to_excel("data/yeast_proteomics_example_scale.xlsx")




# TODO
# test the new coefficient???
coefficient0 = 6.55e9
abundance_ex_corrected["copy_per_cell"] = abundance_ex_corrected["abundance_mmol"] * coefficient0 # Molecular/cell
abundance_ex_corrected.to_excel("data/yeast_proteomics_example_scale.xlsx")

