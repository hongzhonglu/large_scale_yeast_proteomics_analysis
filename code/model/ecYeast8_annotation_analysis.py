# Import packages

import os    ##for directory
import sys
import pprint
from cobra.io import load_matlab_model

os.chdir('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code')
#sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
#pprint.pprint(sys.path)

# import self function
from src.mainFunction import *

# input the subsystem information
gem_dataframe = pd.read_excel('/Users/xluhon/Documents/GitHub/model_correction/result/yeastGEM_with subsystem.xlsx')


# input ecGEM
ecYeast = load_matlab_model('/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat')
#ecYeast = correctSomeWrongFormat(ecYeast)
#produce the dataframe for the metabolites and the rxn
gem_met_nov = produceMetaboliteList(ecYeast)
gem_rxn_nov = produceRxnList(ecYeast)
gem_gene_nov = produceGeneList(ecYeast)



# add subsytem
gem_rxn_nov['subsystem'] = [None]*len(gem_rxn_nov['equation'])
#obtain the transport reaction
gem_rxn_nov['subsystem'] = transport(gem_rxn_nov['formula'].tolist(), gem_rxn_nov['subsystem'].tolist())










