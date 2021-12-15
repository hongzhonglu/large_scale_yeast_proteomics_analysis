# Based on yeast-GEM 8.5

from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model

# input model information
# firstly remove column of "#"
yeast_gem_gene = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="GENES")
yeast_gem = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="RXNS")
# extract rxn and gene
yeast_gem0 = yeast_gem[["ID","GENE ASSOCIATION"]]
yeast_gem0.columns = ["ID", "GPR"]

# read the model
GEM_yeast = read_sbml_model('/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xml')


# run the model
# set the bounds
GEM_yeast.reactions.get_by_id("r_1714").bounds = (-3.255,-3.255) #glu
GEM_yeast.reactions.get_by_id("r_1808").bounds = (0, 0) #glyc
GEM_yeast.reactions.get_by_id("r_2033").bounds = (0.047, 0.047) #pyr
GEM_yeast.reactions.get_by_id("r_2056").bounds = (0,0) #succ
GEM_yeast.reactions.get_by_id("r_1634").bounds = (1.397,1.397) #ac
GEM_yeast.reactions.get_by_id("r_1697").bounds = (7.633,7.633) #co2

# optimize
solution = GEM_yeast.optimize()
GEM_yeast.summary()
fluxes = pd.DataFrame(solution.fluxes)






