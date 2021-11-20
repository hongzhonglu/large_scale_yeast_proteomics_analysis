# this script is to try combine compartment annotation with protein reference abundance and size
# Hongzhong Lu
# 2021-11-20


import os    ##for directory
import numpy as np
import pandas as pd
import math
import sys
import statistics

# os.chdir('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code')
sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
# import self function
from mainFunction import *

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

# input the pro compartment annotation data
pro_location = pd.read_excel("result/gene_compartment_mapping.xlsx")
pro_location = pro_location[['DBID', 'Systematic_name','GO_Name', 'Annot_Type', 'compartment']]

# input the pro abundance from different sources
# should further explore how to scale protein abundance to get absolute concentrations
# input data from SGD
pro_abundance = pd.read_csv("data/sce_protein_abundance_sgd.tsv", sep='\t')
pro_abundance = pro_abundance[pro_abundance["Abundance_median"].notna()]

# with the mitochondrion as an example
genes_mitochondrion = getGeneListFromLocation(pro_location, 'mitochondrion')

# change it as a dataframe
combine_df = pd.DataFrame({"gene": genes_mitochondrion})
combine_df["Volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],combine_df["gene"])
combine_df["section_area"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],combine_df["gene"])
combine_df["abundance"] = singleMapping(pro_abundance["Abundance_median"],pro_abundance['Systematic_name'],combine_df["gene"])

# for the protein without abundance, use the median value from this group.
# calculate the abundance median value
abundance0 = combine_df["abundance"].tolist()
abundance1 = [x for x in abundance0 if np.isnan(x) == False]
abundance_median = statistics.median(abundance1)
abundance_update = []
for x in abundance0:
    if np.isnan(x) == False:
        x0 = x
    else:
        x0 = abundance_median
    abundance_update.append(x0)
combine_df["abundance_update"] = abundance_update

# calculate the size of all proteins in mitochondrion
# 1 纳米(nm)=0.001 微米(um)
total_volume = sum(combine_df["abundance_update"]*combine_df["Volume"])
# change nm^3 into um^3
total_volume_um = total_volume/1e9



