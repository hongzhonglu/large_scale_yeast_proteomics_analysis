# Calcualte the absolute protein abundance from each organelle
# 2022-03-14

import matplotlib.pyplot as plt
from src.protein_process import *

# read the protein abundance files
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx") # all the current datasets


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