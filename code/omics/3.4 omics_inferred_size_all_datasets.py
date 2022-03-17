# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
column20 = list(protein_copy_all1.columns)
column20 = [x for x in column20 if x != "gene"]


# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]


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
result2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")



# only analyze the jianye datasets under C limitation?
all_columns = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in all_columns if "_M" in x]
protein_copy_all1 = protein_copy_all1[Sample_ID_select+["gene"]]
compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
all_compartment = list(compartment.keys())
# compartment curation
# here for the following two membrane, only the compartment annotation with experimental evidence is used.
# also for the fungal-type vacuole membrane, some proteins belong to the metabolic enzymes with high abundance were removed.
# all_compartment = ['plasma membrane']
gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
# all_compartment = ['fungal-type vacuole membrane']
gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
result1 = pd.DataFrame({"compartment": all_compartment})
result2 = pd.DataFrame({"compartment": all_compartment})
for col0 in Sample_ID_select:
    print(col0)
    value1=[]
    value2=[]
    for y in all_compartment:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        if y == "plasma membrane":
            genes_select = gene_plasma_membrane["gene"].tolist()# for the test
        elif y == "fungal-type vacuole membrane":
            genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
            genes_select = [x for x in genes_select if x not in ["YAL005C","YLL024C"]] # remove two genes for fungal type vacuole membrane
        else:
            genes_select = compartment[y]
        pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
        if pro_abundance1 is "no_abundance":
            value1.append(None)
            value2.append(None)
        else:
            x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
            value1.append(x)
            value2.append(S)
    result1[col0] = value1
    result2[col0] = value2
result1.to_excel("data/proteomics/volume_size_across_compartment_jianye_C_limitation.xlsx")
result2.to_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")


















# fungal-type vacuole
#import seaborn as sns
#ss = pro_abundance1
#ss = ss[~ss["molecular/cell"].isna()]
#ss = ss[ss["molecular/cell"]>=60000]
# analyze the density
#sns.displot(ss, x="molecular/cell")








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
protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")
glucose_transporter = pd.read_excel("data/glucose_transporter_abundance_check.xlsx")
glucose_transporter = glucose_transporter[glucose_transporter['Note'].isna()]
omics_glucose_transporter = protein_abundance[protein_abundance['all_gene'].isin(glucose_transporter['gene'])]
omics_glucose_transporter.to_excel("data/omics_glucose_transporter.xlsx")

protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx")
glucose_transporter_copy = protein_copy[protein_copy['gene'].isin(glucose_transporter['gene'])]
glucose_transporter_copy['section_area'] = singleMapping(glucose_transporter['section_area'],glucose_transporter['gene'],glucose_transporter['gene'])
# get the section area
glucose_transporter_copy.to_excel("data/glucose_transporter_copy.xlsx")






