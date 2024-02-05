# -*- coding: utf-8 -*-
'''this code is to read latest yeastGEM, estabolish the gene-protein-reaction relation, then the protein information can be merged into this dataframe
based on the geneID mapping
12th, Nov, 2018
Hongzhong Lu
'''

# Import packages
import pandas as pd
import os    ##for directory
from cobra.io import read_sbml_model
from src.mainFunction import *

# input the subsystem information
gem_dataframe = pd.read_excel('/Users/xluhon/Documents/GitHub/model_correction/result/yeastGEM_with subsystem.xlsx')

# input yeast8 for every update from yeastGEM repo
GEM_nov = read_sbml_model('/Users/xluhon/Documents/GitHub/cobrapy/data/yeastGEM_nov.xml')
#GEM_nov= correctSomeWrongFormat(GEM_nov) # not needed any more for new version of python
#produce the dataframe for the metabolites and the rxn

gem_rxn_nov = produceRxnList(GEM_nov)


#establish rxn-gene mapping
proYeast_DataFrame = getRXNgeneMapping(gem_rxn_nov['rxnID'], gem_rxn_nov['GPR'])
proYeast_DataFrame0 = pd.merge(proYeast_DataFrame, gem_rxn_nov, on='rxnID')




