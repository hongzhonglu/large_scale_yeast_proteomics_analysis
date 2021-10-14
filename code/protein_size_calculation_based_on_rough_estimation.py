# this script is to calculate the protein parameters based on rough estimation
# Hongzhong Lu
# 2021-10-14


import os    ##for directory
import numpy as np
import pandas as pd
import math

# Input the datasets
data_merge = pd.read_csv("data/sce_protein_weight.tsv",sep='\t')

# calculate the size of proteins
radius_list = []
volume_list = []
surface_list = []
for i, x in data_merge.iterrows():
    print(x["proteins_molecular_weight"])
    molecular_weight = x["proteins_molecular_weight"]
    radius = 0.066*molecular_weight**(1/3)
    radius_list.append(radius) # unit is the nm!!
    volume = 4/3*math.pi*radius**3
    volume_list.append(volume)
    surface = 4*math.pi*radius**2
    surface_list.append(surface)


data_merge["radius"] = radius_list
data_merge["volume"] = volume_list
data_merge["surface"] = surface_list

data_merge.to_excel("result/sce_protein_size_rough.xlsx")









