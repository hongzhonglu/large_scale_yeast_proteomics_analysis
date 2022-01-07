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
abundance_ex["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex["genes"])



# Some species process. It could find that some rows with two proteins. Here assume the two proteins are identical in abundance，taking half of total protein abundance.
abudance_ex_1 = abundance_ex[~abundance_ex["MW_Kda"].isna()]
abudance_ex_check = abundance_ex[abundance_ex["MW_Kda"].isna()]
gene0 = []
abundance0 =[]
column0 = abudance_ex_check.columns
abudance_ex_check11 = abudance_ex_check.copy()
gene0=[]
for x in column0:
    if "genes" in x:
        ss = abudance_ex_check[x].tolist()
        for v in ss:
            v1 = v.split("; ")
            gene0= gene0+v1
    else:
        pass

abudance_ex2 = pd.DataFrame({"genes": gene0})


for x in column0:
    if "g/gDW" in x:
        ss = abudance_ex_check[x].tolist()
        abudance0 = []
        for v in ss:
            v1 = [v/2]*2
            abudance0 = abudance0 + v1
        abudance_ex2[x] = abudance0

    else:
        pass





abudance_ex2["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abudance_ex2["genes"])
# combine the above two datasets
abundance_ex_corrected = pd.concat([abudance_ex_1, abudance_ex2], axis=0)
# change the unit from g/g into mmol/g

abundance_ex_corrected["Glucose_phase(mmol/gDW)"] = abundance_ex_corrected["Glucose_phase(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
abundance_ex_corrected["Diauxic_shift(mmol/gDW)"] = abundance_ex_corrected["Diauxic_shift(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
abundance_ex_corrected["Ethanol_phase(mmol/gDW)"] = abundance_ex_corrected["Ethanol_phase(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
new_columns = ['genes',"Glucose_phase(mmol/gDW)","Diauxic_shift(mmol/gDW)","Ethanol_phase(mmol/gDW)"]
abundance_ex_corrected1 = abundance_ex_corrected[new_columns]
#one_strange = abundance_ex_corrected[abundance_ex_corrected['genes']=='YMR142C']
abundance_ex_corrected1.to_excel("data/proteomics/omics_Francesca_scale.xlsx")






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


# input data from paxdb
protein_copy = pd.read_csv("data/proteomics/abundance_table_paxdb.csv")
protein_copy['molecular/cell'] = protein_copy['abundance']*80 # here assume the number of total protein moleculars is 80 million to do the scale!
protein_copy = protein_copy[['gene','molecular/cell']]
protein_copy.columns = ['gene','molecular/cell_paxdb']
protein_copy.to_excel("data/proteomics/abundance_table_paxdb_scale.xlsx", index=False)



# input data from other lab
# this data is from Proteome overabundance enables respiration but limitation onsets carbon overflow
# the unit is molecular/pgDCW, need change it as mmol/gDW
protein_abundance = pd.read_excel("data/proteomics/proteomics_Rahul_2020.xlsx")
colnames = protein_abundance.columns
protein_abundance1 = protein_abundance[colnames[1:]]
colnames0 = colnames[2:]
protein_abundance2 = protein_abundance[colnames[2:]]
protein_abundance3 = protein_abundance[[colnames[1]]]
protein_abundance3.columns = ['gene']
for x in colnames0:
    print(x)
    protein_abundance3[x] = protein_abundance[x]/6.022e23*1000*1e12

protein_abundance3.to_excel("data/proteomics/proteomics_Rahul_2020_scale.xlsx", index=False)









