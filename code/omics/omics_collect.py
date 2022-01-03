import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *


# Part 1 Collect all the data in the unit of mmol/gDW
# input the Jianye's data
growth2 = [0.027, 0.044, 0.102, 0.152, 0.214, 0.254, 0.284, 0.334, 0.379, 0.43]
all_dilution = []
for i in growth2:
    if i < 0.43:
        print(i)
        string0 = "D=" + str(i)
        all_dilution.append(string0)
    else:
        break
# input the measured values
omics_jianye = pd.read_csv("data/proteomics/Omics_from_Jianye.csv")
columns0 = list(omics_jianye.columns)
columns0 = [x for x in columns0 if "RNA" not in x]
columns0 = [x for x in columns0 if "XIA" not in x]
omics_jianye0 = omics_jianye[columns0]
new_columns0 = ['Accession','Gene'] + [x + "_M" for x in all_dilution]
omics_jianye0.columns = new_columns0
omics_jianye1 = omics_jianye0[['Accession','Gene']]

for x in new_columns0:
    if "_M" in x:
        print(x)
        ss1 = omics_jianye0[x]*1e-09
        omics_jianye1[x] = list(ss1)
# id mapping
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
omics_jianye1['gene'] = singleMapping(id_mapping['GeneName'], id_mapping['Entry'], omics_jianye1['Accession'])



# input the tao's data
# input the measured values
omics_tao = pd.read_csv("data/proteomics/Omics_from_tao.csv")
columns0 = list(omics_tao.columns)
columns0 = [x for x in columns0 if "RNA" not in x]
omics_tao0 = omics_tao[columns0]

omics_tao1 = omics_tao0[['Accession','Gene']]

for x in columns0:
    if "prot" in x:
        print(x)
        ss1 = omics_tao0[x]*1e-09
        omics_tao1[x] = list(ss1)
# id mapping
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
omics_tao1['gene'] = singleMapping(id_mapping['GeneName'], id_mapping['Entry'], omics_tao1['Accession'])



# input the Carl's data under max growth
omics_carl = pd.read_excel("data/proteomics/data_PNAS_2021_scale.xlsx")



# input the Francesca's data under max growth
omics_francesca = pd.read_excel("data/proteomics/omics_Francesca_scale.xlsx")



# input the johan's data under max growth
omics_johan = pd.read_excel("data/proteomics/omics_johan.xlsx")





# combine data from different source?
all_gene = set(omics_carl['gene'].tolist()) | set(omics_francesca['genes'].tolist()) | set(omics_tao1['gene'].tolist()) | set(omics_jianye1['gene'].tolist()) | set(omics_johan['gene'].tolist())
all_gene = list(set(all_gene))
new_df = pd.DataFrame({"all_gene": all_gene})
df_combine = pd.merge(left=new_df, right=omics_francesca, left_on=['all_gene'], right_on=['genes'], how="left")
df_combine = df_combine[['all_gene','mmol/gDW']]
df_combine.columns = ['all_gene','mmol/gDW_francesca']

df_combine1 = pd.merge(left=df_combine, right=omics_carl, left_on=['all_gene'], right_on=['gene'], how="left")
df_combine1 = df_combine1[['all_gene','mmol/gDW_francesca', 'mmol/gDW']]
df_combine1.columns = ['all_gene','mmol/gDW_francesca','mmol/gDW_carl']

df_combine2 = pd.merge(left=df_combine1, right=omics_tao1, left_on=['all_gene'], right_on=['gene'], how="left")
df_combine3 = pd.merge(left=df_combine2, right=omics_jianye1, left_on=['all_gene'], right_on=['gene'], how="left")
df_combine4 = pd.merge(left=df_combine3, right=omics_johan, left_on=['all_gene'], right_on=['gene'], how="left")



# get the new column
all0 = df_combine4.columns
new_columns00=[]
for x in all0:
    print(x)
    if "D=" in x:
        new_columns00.append(x)

    elif "prot." in x:
        new_columns00.append(x)

    elif "mmol" in x:
        new_columns00.append(x)
    elif "all_gene" in x:
        new_columns00.append(x)
    else:
        pass

omics_combine = df_combine4[new_columns00]
omics_combine.to_excel("data/proteomics/omics_measured_combine.xlsx")





# Part 2 Collect all the data in the molecular/cell
# Generally, there are three sources, SGD, cell system and another paper. The SGD data use the median value from cell system.
# input data from SGD
# pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')

# input data from cell system, 2018
pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
pro_abundance = pro_abundance[["Systematic Name","Mean molecules per cell","Median molecules per cell"]]
pro_abundance.columns = ["gene", "Mean molecules per cell_cell_system_2018","Median molecules per cell_cell_system_2018"] # protein abundance per cell


# input data from cell reports 2017
pro_abundance2 = pd.read_excel("data/proteomics/protein_copy_cell_report_2017.xlsx")




# combine data from different source?
protein_copy = pd.merge(left=pro_abundance2, right=pro_abundance, left_on=['gene'], right_on=['gene'], how="left")
protein_copy.to_excel("data/proteomics/protein_copy_combine.xlsx")














# application 1 - check the glucose transporter abundance
glucose_transporter = pd.read_excel("data/glucose_transporter_abundance_check.xlsx")
glucose_transporter = glucose_transporter[glucose_transporter['Note'].isna()]
omics_glucose_transporter = omics_combine[omics_combine['all_gene'].isin(glucose_transporter['gene'])]
omics_glucose_transporter.to_excel("data/omics_glucose_transporter.xlsx")


glucose_transporter_copy = protein_copy[protein_copy['gene'].isin(glucose_transporter['gene'])]
glucose_transporter_copy['section_area'] = singleMapping(glucose_transporter['section_area'],glucose_transporter['gene'],glucose_transporter['gene'])
# get the section area
glucose_transporter_copy.to_excel("data/glucose_transporter_copy.xlsx")






