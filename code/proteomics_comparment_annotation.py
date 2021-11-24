# this script is to process compartment datasets
# 2021-11-16

# here the compartment annotation is mainly from SGD and MitoMiner

import os    # for directory
import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt

# function
def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result

# define a function to extract gene list based on its compartment information
def getGeneListFromLocation(gene_location_annotation, location):
    #location = 'mitochondrial envelope'
    gene_subset = gene_location_annotation[gene_location_annotation["GO_Name"]==location]
    gene_list = list(set(gene_subset["Systematic_name"].tolist()))
    return gene_list



# Input the datasets from paxDB
compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')

# extract compartment
compartment.columns = ['DBID', 'Systematic_name', 'Organism','Standard_name', 'Gene_name', 'GO_Qualifier', 'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']

compartment1 = compartment[compartment["GO_Namespace"]=="cellular_component"]
compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
gene_all =list(set(compartment2["DBID"].tolist()))
GO_component_all = list(set(compartment2["GO_Name"].tolist()))
# firstly remove some general cellular component
compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]
# here we need firstly analyze the cellular component information and then filter some unused one?
component_analysis = compartment2['GO_Name'].value_counts()
component_analysis.to_excel("data/all_component_analysis.xlsx")

# re-input the manual check result
component_manual_check = pd.read_excel("data/all_component_manual_check.xlsx")
component_manual_check = component_manual_check[component_manual_check["manual_check"].notna()]

compartment2["compartment"] = singleMapping(component_manual_check["manual_check"],component_manual_check["GO_Name"],compartment2["GO_Name"])

# then all genes could be divided into two groups, with and without compartment
gene_g1 = compartment2[compartment2["compartment"].notna()]


# for the group2, we put it in cytoplasm.
gene_g2 = compartment2[compartment2["compartment"].isna()]
# analyze the result
geneID_g1 = list(set(gene_g1["DBID"].tolist()))
gene_g2_filter1 = gene_g2[~gene_g2["DBID"].isin(geneID_g1)]

# here we assume if we cannot find a compartment, this protein will be in cytoplasm
gene_g2_filter1["compartment"] = "cytoplasm"
len(set(gene_g2_filter1["DBID"].tolist()))


# combine gene_g1 and gene_g2_filter1
gene_combine = pd.concat([gene_g1, gene_g2_filter1], axis=0)


# further analyze the result
gene_combine.to_excel("result/gene_compartment_mapping.xlsx")
#genes_mitochondrion = getGeneListFromLocation(gene_combine, 'mitochondrion')




# Creating a bar chart with the parameters
# first remove the duplications
# this step is not right!! This is because that the column of compartment is from manual check, maybe not right in some protein compartment annotations.

gene_combine["combine"] = gene_combine["DBID"] + "@" + gene_combine["compartment"]
gene_combine1 = gene_combine.drop_duplicates(subset='combine', keep="last")
component_analysis2 = gene_combine1['compartment'].value_counts()
component_analysis2.to_excel("result/main_compartment_statistical_analysis.xlsx")

# plot
x_value = list(component_analysis2.index)
y_value = list(component_analysis2.values)
plt.figure(figsize=(4,3))
plt.bar(x_value, y_value, width=0.7, bottom=50, align='edge')
plt.title('')
plt.xlabel('Compartment', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(rotation=90)
plt.show()
#plt.savefig("result/gene_compartment_analysis.pdf")


# special analysis of mitochrondria
# for the group1, we firstly look into mitochrondria
gene_g1_m = gene_g1[gene_g1["compartment"]=="mitochondrion/mitochondrial membrane"]
gene_g1_others = gene_g1[~(gene_g1["compartment"]=="mitochondrion/mitochondrial membrane")]
# try to remove duplicates in mitochrondria
gene_g1_m.to_excel("data/gene_g1_m.xlsx")
genes_mitochondrion = getGeneListFromLocation(gene_g1_m, 'mitochondrion')
genes_mito_om = getGeneListFromLocation(gene_g1_m, 'mitochondrial outer membrane')
genes_mito_im = getGeneListFromLocation(gene_g1_m, 'mitochondrial inner membrane')
genes_mito_is = getGeneListFromLocation(gene_g1_m, 'mitochondrial intermembrane space')
genes_mito_matrix = getGeneListFromLocation(gene_g1_m, 'mitochondrial matrix')
genes_mito_membrane = getGeneListFromLocation(gene_g1_m, 'mitochondrial membrane')
genes_mito_envelope = getGeneListFromLocation(gene_g1_m, 'mitochondrial envelope')
# it is shown some genes in the envelope, but not belong to om or im
# only three proteins were not assigned as om or im,
# so we only assume they are in mitochondrion.
# It is also shown that some proteins were annotated as membrane but not at im, om and is.
# gene_need_check = list(set(genes_mito_envelope)-set(genes_mito_om + genes_mito_im))
# gene_need_check1 = list(set(genes_mito_membrane) - set(genes_mito_om + genes_mito_im + genes_mito_is))
genes_cytoplasm = getGeneListFromLocation(compartment2, 'cytoplasm')
genes_cytosol = getGeneListFromLocation(compartment2, 'cytosol')