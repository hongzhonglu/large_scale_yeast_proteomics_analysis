# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.protein_process import *



# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]


# compartment
compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
all_compartment = list(compartment.keys())


# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
Sample_ID_select = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]


# recalculation based on the curated protein abundances
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
result1.to_excel("data/proteomics/volume_size_across_compartment.xlsx")
result2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")



# go term
go_term = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")
all_go_term = list(go_term.keys())
result1 = pd.DataFrame({"go_term": all_go_term})
result2 = pd.DataFrame({"go_term": all_go_term})

for col0 in Sample_ID_select:
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
            x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
            value1.append(x)
            value2.append(S)
    result1[col0] = value1
    result2[col0] = value2

result1.to_excel("data/proteomics/volume_size_across_go_term.xlsx")
result2.to_excel("data/proteomics/membrance_size_across_go_term.xlsx")


