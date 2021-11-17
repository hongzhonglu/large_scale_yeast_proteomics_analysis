# this script is to process compartment datasets
# 2021-11-16

# here the comparment annotation is mainly from SGD and MitoMiner

import os    # for directory
import numpy as np
import pandas as pd
import math

# Input the datasets from paxDB
compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')

# extract compartment
compartment.columns = ['DBID', 'Systematic_name', 'Organism','Standard_name', 'Gene_name', 'GO_Qualifier', 'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']

compartment1 = compartment[compartment["GO_Namespace"]=="cellular_component"]
compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
gene_all =list(set(compartment2["DBID"].tolist()))
GO_component_all = list(set(compartment2["GO_Name"].tolist()))

# here we need firstly analyze the cellular component information and then filter some unused one?
component_analysis = compartment2['GO_Name'].value_counts()



compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]

component_analysis = compartment2['GO_Name'].value_counts()

component_analysis.to_excel("data/all_component_analysis.xlsx")




