# Based on yeast-GEM 8.5

from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from src.model_process import *
from src.mainFunction import *

# input model information
# firstly remove column of "#"
yeast_gem = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", usecols="A:P", sheet_name="RXNS")


# check which reaction has no subsystems
subsystem = pd.read_csv('data/subsystem/Rxn_unique_subsystem.tsv', sep='\t')
new_rxn = yeast_gem[~yeast_gem['ID'].isin(subsystem['Abbreviation'])]
old_rxn = yeast_gem[yeast_gem['ID'].isin(subsystem['Abbreviation'])]
old_rxn['subsystem_unique'] = singleMapping(subsystem['Subsystem_new'], subsystem['Abbreviation'], old_rxn['ID'])



# extract the subsytem for these new rxns
new_rxn["subsystem_unique"] = [None]*new_rxn.shape[0]
new_subsystem = new_rxn["subsystem_unique"].tolist()
s0 = exchange(s1=new_rxn['EQUATION'], subystem=new_subsystem)
s1 = transport(s1=new_rxn['EQUATION'], subsysem=s0)
new_rxn["subsystem_unique"] = s1
new_rxn["subsystem_unique"][new_rxn["ID"]=="r_4598"] = "growth"
new_rxn["subsystem_unique"][new_rxn["ID"]=="r_4599"] = "growth"


for i, x in new_rxn.iterrows():
    if "FA ester" in x["NOTE"] and x['subsystem_unique'] is None:
        new_rxn["subsystem_unique"][i] = 'FA ester pathway'
    else:
        new_rxn["subsystem_unique"][i] = new_rxn["subsystem_unique"][i]


# combine the subsystem
all_rxn_with_subsystem = pd.concat([old_rxn,new_rxn],axis=0)


# unify the format of transporter reaction:
all_rxn_with_subsystem['subsystem_unique'] = transport(s1=all_rxn_with_subsystem['EQUATION'], subsysem=all_rxn_with_subsystem['subsystem_unique'])
all_rxn_with_subsystem.to_csv("data/subsystem/Rxn_unique_subsystem_v2.tsv", sep="\t")








# read the model
GEM_yeast = read_sbml_model('/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xml')
gene_GEM = getALLGEMgene()












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






