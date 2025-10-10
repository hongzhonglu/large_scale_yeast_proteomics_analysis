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







############## datasets preprocess ###################################
# This original dataset is from www.pnas.org/cgi/doi/10.1073/pnas.1918216117
# Note: the full dataset is used
abundance_ex = pd.read_excel("data/proteomics/omics_Francesca.xlsx")
# remove one outlier data point
# abundance_ex['Glucose_phase_rep1(g/gDW)'][abundance_ex["gene"]=="YMR142C"] = None
abundance_ex = abundance_ex[abundance_ex["gene"] != "YMR142C"]

abundance_ex["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex["gene"])
# Some species process. It could find that some rows with two proteins. Here assume the two proteins are identical in abundance，taking half of total protein abundance.
abudance_ex_1 = abundance_ex[~abundance_ex["MW_Kda"].isna()]
abudance_ex_check = abundance_ex[abundance_ex["MW_Kda"].isna()]
gene0 = []
abundance0 =[]
column0 = abudance_ex_check.columns
abudance_ex_check11 = abudance_ex_check.copy()
gene0=[]
for x in column0:
    if "gene" in x:
        ss = abudance_ex_check[x].tolist()
        for v in ss:
            v1 = v.split("; ")
            gene0= gene0+v1
    else:
        pass

abudance_ex2 = pd.DataFrame({"gene": gene0})


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


abudance_ex2["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abudance_ex2["gene"])
# combine the above two datasets
abundance_ex_corrected = pd.concat([abudance_ex_1, abudance_ex2], axis=0)
# change the unit from g/g into mmol/g
"""
abundance_ex_corrected["Glucose_phase(mmol/gDW)"] = abundance_ex_corrected["Glucose_phase(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
abundance_ex_corrected["Diauxic_shift(mmol/gDW)"] = abundance_ex_corrected["Diauxic_shift(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
abundance_ex_corrected["Ethanol_phase(mmol/gDW)"] = abundance_ex_corrected["Ethanol_phase(g/gDW)"]/abundance_ex_corrected["MW_Kda"]# #mmol/g biomass
"""
column_ss = list(abundance_ex.columns)[1:10]
for xx in column_ss:
    abundance_ex_corrected[xx] = abundance_ex_corrected[xx] / abundance_ex_corrected["MW_Kda"]  # #mmol/g biomass

new_columns = ['gene'] + column_ss
abundance_ex_corrected1 = abundance_ex_corrected[new_columns]
#one_strange = abundance_ex_corrected[abundance_ex_corrected['gene']=='YMR142C']
abundance_ex_corrected1.to_excel("data/proteomics/omics_Francesca_scale.xlsx")





############## datasets preprocess ###################################
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




############## datasets preprocess ###################################
# input data from paxdb
protein_copy = pd.read_csv("data/proteomics/abundance_table_paxdb.csv")
protein_copy['molecular/cell'] = protein_copy['abundance']*80 # here assume the number of total protein moleculars is 80 million to do the scale!
protein_copy = protein_copy[['gene','molecular/cell']]
protein_copy.columns = ['gene','molecular/cell_paxdb']
protein_copy.to_excel("data/proteomics/abundance_table_paxdb_scale.xlsx", index=False)





############## datasets preprocess ###################################
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






############## datasets preprocess ###################################
# input the Jianye's data
# in the original Jianye datasets, the unit is protein copies/cell, so no additional unit conversion is needed!
# Note: in the current calculation by Jianye, the assume cell mass is 13 pg.
# Note: the part of data has been updated
omics_jianye_original = pd.read_excel("data/proteomics/proteomics_Jianye_original.xlsx")
# first update isoform of proteins and change it as a single protein
ID_new = []
for i, x in omics_jianye_original.iterrows():
    #print(i)
    ss = x["Majority protein IDs"]
    if "-2" in x["Majority protein IDs"]:
        print(x["Majority protein IDs"])
        ss1 = ss.replace("-2","").replace("-3","")
        ss2 = ss1.split(";")[0]
        print(ss2)
        ID_new.append(ss2)
    else:
        ID_new.append(ss)

omics_jianye_original1 = omics_jianye_original
omics_jianye_original1["Majority protein IDs"] = ID_new


abudance_jianye_1 = omics_jianye_original1[~omics_jianye_original1["Majority protein IDs"].str.contains(";")]
abudance_jianye_check = omics_jianye_original1[omics_jianye_original1["Majority protein IDs"].str.contains(";")]
gene0 = []
abundance0 =[]
column0 = abudance_jianye_check.columns
abudance_jianye_check11 = abudance_jianye_check.copy()
gene0 = []
num0 = []
for x in column0:
    if "Majority protein IDs" in x:
        ss = abudance_jianye_check[x].tolist()
        for v in ss:
            v1 = v.split(";")
            gene0= gene0+v1
            num0.append(len(v1))
    else:
        pass

abudance_jianye2 = pd.DataFrame({"Majority protein IDs": gene0})


for x in column0:
    if "S" in x:
        print(x)
        ss = abudance_jianye_check[x].tolist()
        abudance0 = []
        for v,n in zip(ss,num0):
            v1 = [v/n]*n
            abudance0 = abudance0 + v1
        abudance_jianye2[x] = abudance0

    else:
        pass


column1 = [x for x in column0 if x !="Gene Name"]
abudance_jianye_1=abudance_jianye_1[column1]
abundance_jianye_corrected = pd.concat([abudance_jianye_1, abudance_jianye2], axis=0)
#abundance_jianye_corrected['gene'] = singleMapping(id_mapping['GeneName'], id_mapping['Entry'], abundance_jianye_corrected['Majority protein IDs'])
abundance_jianye_corrected['gene'] = multiMapping(id_mapping['GeneName'], id_mapping['Entry'], abundance_jianye_corrected['Majority protein IDs'])

column2 = ["gene"] + column1[1:]
abundance_jianye_corrected = abundance_jianye_corrected[column2]

# from multiMapping function, it could find one uniprot ID could have multiple locus gene ID
column_jianye2 = list(abundance_jianye_corrected.columns)
column_jianye2 = column_jianye2[1:28]
pd_null = pd.DataFrame()
for x in column_jianye2:
    select0 = ["gene", x]
    select_value = abundance_jianye_corrected[select0]
    select_value1 = splitAbundance(select_value)
    pd_null = pd.concat([pd_null, select_value1], axis=1)

_, i = np.unique(pd_null.columns, return_index=True)
omics_jianye_2 = pd_null.iloc[:, i]

omics_jianye_2.to_excel("data/proteomics/protein_copy_jianye.xlsx", index=False)

# change it as abundance, from protein copy/cell to mmol/gDW
coefficient1 = 7.8298e9
abundance_jianye3 = omics_jianye_2.copy()
for x in column_jianye2:
    print(x)
    if x != "gene":
        abundance_jianye3[x] = omics_jianye_2[x]/coefficient1
    else:
        continue
abundance_jianye3 = abundance_jianye3[["gene"]+column_jianye2]
abundance_jianye3.to_excel("data/proteomics/abundance_jianye.xlsx", index=False)

############## datasets preprocess ###################################
# this is dataset from johan
# assume the 0.48 g protein per gram of biomass
# the original unit is (g protein/g proteome), change it as mmol protein/gDW
# In the last version, the proteomics is from johan
abundance_johan = pd.read_excel("data/proteomics/proteomics_johan_original.xlsx")
abundance_johan["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_johan["gene"])
column_johan = list(abundance_johan.columns)[2:12]
for x in column_johan:
    print(x)
    if x != "gene":
        abundance_johan[x] = abundance_johan[x]/abundance_johan["MW_Kda"]*0.48
    else:
        continue

abundance_johan = abundance_johan[["gene"]+column_johan]
abundance_johan.to_excel("data/proteomics/omics_johan_all_covered.xlsx", index=False)


############## datasets preprocess ###################################
# this is dataset from Kate
# the unit is protein molecules per pg CDW
abundance_kate = pd.read_excel("data/proteomics/datasets_kate_2020.xlsx")
sample_info_kate = pd.read_excel("data/proteomics/datasets_kate_2020.xlsx", sheet_name="sample_information")
sample_info_kate = sample_info_kate.transpose()
sample_info_kate1 = sample_info_kate.iloc[1: , :]
sample_info_kate1.columns = sample_info_kate.iloc[0]
sample_info_kate1["sample_ID"] = list(sample_info_kate1.index)

# change the unit as mmol/gDW
abundance_kate2 = abundance_kate
column_k = abundance_kate.columns
for x in column_k:
    print(x)
    if x != "gene":
        abundance_kate2[x] = abundance_kate[x]/6.022e23*1000*1e12
    else:
        continue

abundance_kate2.to_excel("data/proteomics/abundance_kate.xlsx", index=False)


# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
abundance_kate3 = abundance_kate2.copy()
for x in column_k:
    print(x)
    if x != "gene":
        abundance_kate3[x] = abundance_kate2[x]*coefficient1
    else:
        continue

abundance_kate3.to_excel("data/proteomics/protein_copy_kate.xlsx", index=False)

# statistical analysis of all proteomics datasets
result_df = AllProteomicsAnalysis(pro_df=abundance_kate3)
result_df0 = result_df.transpose()
result_df0["total_copy"] = result_df0["count"]*result_df0["mean"]





############## datasets preprocess ###################################
# input the tao's data
# Absolute protein and mRNA abundances (fmol/mgDW) by rosemery
# input the measured values
omics_tao = pd.read_csv("data/proteomics/Omics_from_tao.csv")
columns0 = list(omics_tao.columns)
columns0 = [x for x in columns0 if "RNA" not in x]
omics_tao0 = omics_tao[columns0]

omics_tao1 = omics_tao0[['Accession','Gene']]

for x in columns0:
    if "prot" in x:
        print(x)
        ss1 = omics_tao0[x]*1e-09    # change the unit as mmol/gDW!
        omics_tao1[x] = list(ss1)
# id mapping
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
omics_tao1['gene'] = multiMapping(id_mapping['GeneName'], id_mapping['Entry'], omics_tao1['Accession'])


# from multiMapping function, it could find one uniprot ID could have multiple locus gene ID
column_tao2 = list(omics_tao1.columns)
column_tao2 = column_tao2[3:45]
pd_null = pd.DataFrame()
for x in column_tao2:
    select0 = ["gene", x]
    select_value = omics_tao1[select0]
    select_value1 = splitAbundance(select_value)
    pd_null = pd.concat([pd_null, select_value1], axis=1)

_, i = np.unique(pd_null.columns, return_index=True)
omics_tao2 = pd_null.iloc[:, i]
omics_tao2 = omics_tao2[["gene"]+column_tao2]
omics_tao2.to_excel("data/proteomics/Omics_from_tao_scale.xlsx", index=False)



############## datasets preprocess ###################################
# input the tao's data published in Nature communication, 2020.
# Absolute protein and mRNA abundances (fmol/mgDW) by rosemery
# input the measured values
omics_tao_nc = pd.read_excel("data/proteomics/Proteomics_NC_2020_rosmary.xlsx")
columns0 = list(omics_tao_nc.columns)
columns0 = [x for x in columns0 if "RNA" not in x]
omics_tao_nc0 = omics_tao_nc[columns0]

omics_tao_nc1 = omics_tao_nc0[['Accession','Gene']]

for x in columns0:
    if "prot" in x:
        print(x)
        ss1 = omics_tao_nc0[x]*1e-09
        omics_tao_nc1[x] = list(ss1)
# id mapping
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
omics_tao_nc1['gene'] = multiMapping(id_mapping['GeneName'], id_mapping['Entry'], omics_tao_nc1['Accession'])


# from multiMapping function, it could find one uniprot ID could have multiple locus gene ID
column_tao2 = list(omics_tao_nc1.columns)
column_tao2 = column_tao2[2:10]
pd_null = pd.DataFrame()
for x in column_tao2:
    select0 = ["gene", x]
    select_value = omics_tao_nc1[select0]
    select_value1 = splitAbundance(select_value)
    pd_null = pd.concat([pd_null, select_value1], axis=1)

_, i = np.unique(pd_null.columns, return_index=True)
omics_tao_nc2 = pd_null.iloc[:, i]
omics_tao_nc2.to_excel("data/proteomics/Omics_from_tao_nc_scale.xlsx", index=False)


# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
omics_tao_nc3 = omics_tao_nc2.copy()
for x in column_tao2:
    print(x)
    if x != "gene":
        omics_tao_nc3[x] = omics_tao_nc3[x]*coefficient1
    else:
        continue

omics_tao_nc3.to_excel("data/proteomics/protein_copy_tao_nc.xlsx", index=False)

# statistical analysis of all proteomics datasets
result_df = AllProteomicsAnalysis(pro_df=omics_tao_nc3)
result_df0 = result_df.transpose()
result_df0["total_copy"] = result_df0["count"]*result_df0["mean"]





############## datasets preprocess ###################################
# this data is sysbio, cell systems, 2017
# the unit is molecular/pgDCW, need change it as mmol/gDW
protein_abundance = pd.read_excel("data/proteomics/omics_from_cell_systems_2017.xlsx")
protein_abundance = protein_abundance.drop_duplicates("gene", keep='first')

colnames = protein_abundance.columns
protein_abundance1 = protein_abundance[colnames[1:]]
colnames0 = colnames[1:]
protein_abundance3 = protein_abundance[[colnames[0]]]
protein_abundance3.columns = ['gene']
for x in colnames0:
    print(x)
    protein_abundance3[x] = protein_abundance[x]/6.022e23*1000*1e12

protein_abundance3.to_excel("data/proteomics/omics_from_cell_systems_2017_scale.xlsx", index=False)

# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
protein_abundance4 = protein_abundance3.copy()
for x in colnames0:
    print(x)
    if x != "gene":
        protein_abundance4[x] = protein_abundance3[x]*coefficient1
    else:
        continue

protein_abundance4.to_excel("data/proteomics/protein_copy_from_cell_systems_2017_scale.xlsx", index=False)
# statistical analysis of all proteomics datasets
result_df = AllProteomicsAnalysis(pro_df=protein_abundance4)
result_df0 = result_df.transpose()
result_df0["total_copy"] = result_df0["count"]*result_df0["mean"]



############## datasets preprocess ###################################
# TODO
# different substrates
# the unit is mmol/gDW
carbon_source1 = pd.read_excel("data/proteomics/Proteome_ref_carbon_source.xlsx", sheet_name="Sce_carbon1")
carbon_source2 = pd.read_excel("data/proteomics/Proteome_ref_carbon_source.xlsx", sheet_name="Sce_carbon2")
colname_s1 = carbon_source1.columns
colname_s1 = [x for x in colname_s1 if "mmol" in x]

colname_s2 = carbon_source2.columns
colname_s2 = [x for x in colname_s2 if "mmol" in x]

carbon_source10 = carbon_source1[['gene'] + colname_s1]
carbon_source20 = carbon_source2[['gene'] + colname_s2]
carbon_source_combine = pd.merge(left=carbon_source10, right=carbon_source20, left_on=['gene'], right_on=['gene'], how="left")
carbon_source_combine = carbon_source_combine.drop(columns="ref_glc_mm_rich_aerobic(mmol/gDW)_y")
carbon_source_combine.pop('ref_glc_mm_rich_aerobic(mmol/gDW)_x')

carbon_source_combine.to_excel("data/proteomics/omics_from_carbon_source_scale.xlsx", index=False)


colnames0 = carbon_source_combine.columns
# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
carbon_source_combine2 = carbon_source_combine.copy()
for x in colnames0:
    print(x)
    if x != "gene":
        carbon_source_combine2[x] = carbon_source_combine[x]*coefficient1
    else:
        continue

carbon_source_combine2.to_excel("data/proteomics/protein_copy_from_carbon_source_scale.xlsx", index=False)

# statistical analysis of all proteomics datasets
result_df = AllProteomicsAnalysis(pro_df=carbon_source_combine2)
result_df0 = result_df.transpose()
result_df0["total_copy"] = result_df0["count"]*result_df0["mean"]
