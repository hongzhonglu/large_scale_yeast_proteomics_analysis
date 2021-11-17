# this script is to process proteomics datasets
# 2021-11-16


import os    ##for directory
import numpy as np
import pandas as pd
import math

# Input the datasets from paxDB
abundance = pd.read_csv("data/info_yeast/abundance_table.csv")

# Get the molecular weight data
mw = pd.read_csv("data/info_yeast/all_peptides_mw.csv")
mw["MW_Kda"] = mw["MW"]/1000

# combine the above data
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

abundance["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance["genes"])

# scale the dataset
coefficient_protein_in_biomass = 0.45
abundance["abundance_w1"] = abundance["MW_Kda"]*abundance["abundance"]
w_total = sum(abundance["abundance_w1"])
abundance["abundance_scale"] = abundance["abundance_w1"]/w_total
abundance["abundance_per_biomass"] = abundance["abundance_scale"]*coefficient_protein_in_biomass # g/g biomass
abundance["abundance_per_biomass2"] = abundance["abundance_per_biomass"]/abundance["MW_Kda"]# #mmol/g biomass

# then how to transform the absolute abundance into one cell?
# according to this paper https://www.pnas.org/content/107/3/999, we can get the density of yeast cell
# assume
density = 1.1029 #1.1029 ± 0.0026 g/mL
# one weight of a cell
Vcell = 82 #µm3
Mcell = 1.1029*82/1e12 # g
Ncell = 1/Mcell #1.1e7

# set a coefficent to transform the mmol/g.Biomass into molecular/cell
# this step should be careful in the late calculation
coefficient0 = 6.022e20/Ncell
abundance["abundance_per_biomass3"] = abundance["abundance_per_biomass2"]*coefficient0 # Molecular/cell













