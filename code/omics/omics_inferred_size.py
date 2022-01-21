# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *


# protein abundance combine and analysis
protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")
protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx")
# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
column0 = protein_abundance.columns
column1 = column0[1:]
protein_abundance1 = protein_abundance[column1]
protein_abundance2 = protein_abundance1.copy()
for x in column1:
    print(x)
    protein_abundance2[x] = protein_abundance1[x]*coefficient1

protein_abundance2['gene'] = protein_abundance['all_gene']
# combine the two dataframe
protein_copy_all = pd.merge(left=protein_abundance2, right=protein_copy, left_on=['gene'], right_on=['gene'], how='outer')
column2 = protein_copy_all.columns
column20 = [x for x in column2 if x is not 'gene']
column21 = ['gene'] + column20
protein_copy_all1 = protein_copy_all[column21]
protein_copy_all1.to_excel("data/proteomics/all_protein_copy.xlsx")
# statistical analysis of all proteomics datasets
AllProteomicsAnalysis(pro_df=protein_copy_all1)







# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]
# input the pro compartment annotation data
# the protein location here is based on manual check, but it contains annotation from computation!!
pro_location = pd.read_excel("result/gene_compartment_mapping.xlsx")
pro_location = pro_location[['DBID', 'Systematic_name','GO_Name', 'Annot_Type', 'compartment']]








# compartment
compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
all_compartment = list(compartment.keys())
result1 = pd.DataFrame({"compartment": all_compartment})
result2 = pd.DataFrame({"compartment": all_compartment})
for col0 in column20:
    print(col0)
    value1=[]
    value2=[]
    for y in all_compartment:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        genes_select = compartment[y]
        pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
        if pro_abundance1 is "no_abundance":
            value1.append(None)
            value2.append(None)
        else:
            x, S = getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)
            value1.append(x)
            value2.append(S)
    result1[col0] = value1
    result2[col0] = value2
result1.to_excel("data/proteomics/volume_size_across_compartment.xlsx")
result2.to_excel("data/proteomics/membrance_size_across_compartment.xlsx")




# change the size calculation as a function, so that no matter which kind of data as input, we can
# calcualte the result
compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
location0 = "mitochondrial inner membrane"
genes_select = compartment[location0]
pro_df = protein_copy_all1
# abundance check
pro_df_test = pro_df[pro_df["gene"].isin(genes_select)]
pro_df_test.to_excel("data/abundance_check_test.xlsx")



# read the model
GEM_yeast = read_sbml_model('/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xml')
gene_GEM = getALLGEMgene()
genes_select_M = list(set(genes_select) & set(gene_GEM))




protein_copy_all2 = pro_df.drop(columns=['gene'])
column20 = protein_copy_all2.columns
volume0 = []
area0 = []
for col0 in column20:
    print(col0)
    pro_abundance = protein_copy_all1[['gene', col0]]
    pro_abundance.columns = ['gene', 'molecular/cell']
    pro_abundance1 = getProAundance(genes_select0=genes_select_M, pro_abundance0=pro_abundance) # update the gene list
    if pro_abundance1 is "no_abundance":
        volume0.append(None)
        area0.append(None)
    else:
        x, S = getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)
        volume0.append(x)
        area0.append(S)
result_out = pd.DataFrame({"sampelID": column20})
result_out["volume"] = volume0
result_out["section_area"] = area0
result_out.to_excel("data/" + location0 + "_proteomics_size.xlsx")















# go term
go_term = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")
all_go_term = list(go_term.keys())
result1 = pd.DataFrame({"go_term": all_go_term})
result2 = pd.DataFrame({"go_term": all_go_term})

for col0 in column20:
    print(col0)
    value1=[]
    value2=[]
    for y in all_go_term:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        genes_select = go_term[y]
        pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
        if pro_abundance1 is "no_abundance":
            value1.append(None)
            value2.append(None)
        else:
            x, S = getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)
            value1.append(x)
            value2.append(S)
    result1[col0] = value1
    result2[col0] = value2

result1.to_excel("data/proteomics/volume_size_across_go_term.xlsx")
result2.to_excel("data/proteomics/membrance_size_across_go_term.xlsx")















# other for specific reactions - check the glucose transporter abundance
glucose_transporter = pd.read_excel("data/glucose_transporter_abundance_check.xlsx")
glucose_transporter = glucose_transporter[glucose_transporter['Note'].isna()]
omics_glucose_transporter = protein_abundance[protein_abundance['all_gene'].isin(glucose_transporter['gene'])]
omics_glucose_transporter.to_excel("data/omics_glucose_transporter.xlsx")


glucose_transporter_copy = protein_copy[protein_copy['gene'].isin(glucose_transporter['gene'])]
glucose_transporter_copy['section_area'] = singleMapping(glucose_transporter['section_area'],glucose_transporter['gene'],glucose_transporter['gene'])
# get the section area
glucose_transporter_copy.to_excel("data/glucose_transporter_copy.xlsx")






