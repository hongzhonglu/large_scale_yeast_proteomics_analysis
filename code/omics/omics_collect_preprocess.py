# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

import sys

# import self function
from src.mainFunction import *
from src.protein_process import *


# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000


# process the experimentally proteomic datasets
# This original dataset is from www.pnas.org/cgi/doi/10.1073/pnas.1918216117
abundance_ex = pd.read_excel("data/proteomics/omics_Francesca.xlsx")
abundance_ex = abundance_ex[["GeneNameOrdered","average(g/gDW)"]]
abundance_ex.columns = ["genes","g/gDW"]
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
abundance_ex_corrected.to_excel("data/proteomics/omics_Francesca_scale.xlsx")






# input latest dataset of Carl
coefficient1 = 6.5789e9
abundance_ex = pd.read_excel("data/proteomics/data_PNAS_2021.xlsx")
abundance_ex['g/gDW'] =(abundance_ex['replicate 1 (g gDW-1)']+ abundance_ex['replicate 2 (g gDW-1)']+ abundance_ex['replicate 3 (g gDW-1)'])/3
abundance_ex=abundance_ex[['Symbol','g/gDW']]
abundance_ex.columns = ['gene','g/gDW']
abundance_ex1 = splitAbundance(pro_df=abundance_ex)
abundance_ex1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex1["gene"])
abundance_ex_check = abundance_ex1[abundance_ex1["MW_Kda"].isna()]
abundance_ex1=abundance_ex1[~abundance_ex1["MW_Kda"].isna()]
abundance_ex1["mmol/gDW"] = abundance_ex1["g/gDW"]/abundance_ex1["MW_Kda"]# #mmol/g biomass
abundance_ex1["molecular/cell"] = abundance_ex1["mmol/gDW"]*coefficient1 # unit transformation
abundance_ex1.to_excel("data/proteomics/data_PNAS_2021_scale.xlsx")

