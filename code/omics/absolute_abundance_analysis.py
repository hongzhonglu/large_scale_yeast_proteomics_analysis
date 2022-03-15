# Calcualte the absolute protein abundance from each organelle
# 2022-03-14

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



# test: check the calculation of abundance sum for gene from one compartment
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all1 = protein_copy_all1[Sample_ID_select+["gene"]]


compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
all_compartment = list(compartment.keys())
result1 = pd.DataFrame({"compartment": all_compartment})
for col0 in Sample_ID_select:
    print(col0)
    value1=[]
    for y in all_compartment:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        genes_select = compartment[y]
        pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
        if pro_abundance1 is "no_abundance":
            value1.append(None)
        else:
            x = getOrganelleAbundance(abundance0=pro_abundance1)
            value1.append(x)
    result1[col0] = value1

result1.to_excel("data/proteomics/protein_abundance_across_compartment.xlsx")
