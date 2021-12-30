# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

import sys
# sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
# import self function
from src.mainFunction import *
from src.protein_process import *


# coeffcient change
# this data maybe not right!
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



# To calculate the structure constraint , we need change the unit g protein per biomass into  molecular per cell !

# process the experimentally proteomic datasets
# This original dataset is from www.pnas.org/cgi/doi/10.1073/pnas.1918216117
abundance_ex = pd.read_excel("data/proteomics/omics_Francesca.xlsx")
abundance_ex = abundance_ex[["GeneNameOrdered","average(g/gDW)"]]
abundance_ex.columns = ["genes","g/gDW"]

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
    v = [x["g/gDW"]/2]*2
    abundance0 = abundance0 + v
abudance_ex_2 = pd.DataFrame({"genes":gene0,"g/gDW":abundance0})
abudance_ex_2["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abudance_ex_2["genes"])

# combine the above two datasets
abundance_ex_corrected = pd.concat([abudance_ex_1, abudance_ex_2], axis=0)
# change the unit from g/g into mmol/g
abundance_ex_corrected["mmol/gDW"] = abundance_ex_corrected["g/gDW"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
one_strange = abundance_ex_corrected[abundance_ex_corrected['genes']=='YMR142C']
abundance_ex_corrected["molecular/cell"] = abundance_ex_corrected["mmol/gDW"] * coefficient0 # Molecular/cell
abundance_ex_corrected.to_excel("data/proteomics/omics_Francesca_scale.xlsx")


# TODO
# test the new coefficient???
coefficient0 = 6.5789e9 # this coefficent is from proteomics quality check to transfer mmol/gDW into molecular/cell
abundance_ex_corrected["molecular/cell"] = abundance_ex_corrected["mmol/gDW"] * coefficient0 # Molecular/cell
abundance_ex_corrected.to_excel("data/proteomics/omics_Francesca_scale.xlsx")






# input latest dataset of Carl
abundance_ex = pd.read_excel("data/proteomics/data_PNAS_2021.xlsx")
abundance_ex['g/gDW'] =(abundance_ex['replicate 1 (g gDW-1)']+ abundance_ex['replicate 2 (g gDW-1)']+ abundance_ex['replicate 3 (g gDW-1)'])/3
abundance_ex=abundance_ex[['Symbol','g/gDW']]
abundance_ex.columns = ['gene','g/gDW']

abundance_ex1 = splitAbundance(pro_df=abundance_ex)
sum(abundance_ex1['g/gDW'])


mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000

abundance_ex1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex1["gene"])

abundance_ex_check = abundance_ex1[abundance_ex1["MW_Kda"].isna()]
abundance_ex1=abundance_ex1[~abundance_ex1["MW_Kda"].isna()]


abundance_ex1["mmol/gDW"] = abundance_ex1["g/gDW"]/abundance_ex1["MW_Kda"]# #mmol/g biomass

abundance_ex1["molecular/cell"] = abundance_ex1["mmol/gDW"] * coefficient0 # Molecular/cell
abundance_ex1.to_excel("data/proteomics/data_PNAS_2021_scale.xlsx")
sum(abundance_ex1["molecular/cell"])
