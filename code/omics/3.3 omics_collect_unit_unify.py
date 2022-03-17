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
protein_copy_all1.to_excel("data/proteomics/all_protein_copy.xlsx", index=False)



