# -*- coding: utf-8 -*-
'''this code is to carry out the model validation analysis based on yeastGEM with different sources data;
new gene and new reaction analysis for Yeast8 compared with Yeast7
5th, Nov, 2018'''

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
ecYeast = correctSomeWrongFormat(ecYeast)
#produce the dataframe for the metabolites and the rxn
gem_met_nov = produceMetaboliteList(ecYeast)
gem_rxn_nov = produceRxnList(ecYeast)
gem_gene_nov = produceGeneList(ecYeast)

# add subsytem
gem_rxn_nov['subsystem'] = [None]*len(gem_rxn_nov['equation'])
#obtain the transport reaction
gem_rxn_nov['subsystem'] = transport(gem_rxn_nov['formula'].tolist(), gem_rxn_nov['subsystem'].tolist())






















#the following is not used!
#remove the transport reaction
gem_rxn_transport = gem_rxn_nov[gem_rxn_nov['subsystem'].str.contains('Transport', na=False)]
gem_rxn_metabolic = gem_rxn_nov[gem_rxn_nov['subsystem'].str.contains('Transport', na=False)==False]
#remove the exchange reaction
gem_rxn_metabolic1 = gem_rxn_metabolic[gem_rxn_metabolic['name'].str.contains('exchange', na=False)==False]
#remove the arm reaction
gem_rxn_metabolic2 = gem_rxn_metabolic1[gem_rxn_metabolic1['rxnID'].str.contains('arm', na=False)==False]
#all metabolic reaction
gpr = gem_rxn_metabolic2['GPR'].tolist()
gpr[1]
#obtain number of metabolic reaction with gpr
rxn_gpr = [x for x in gpr if x !='']
#obtain the number of metabolic reaction without gpr
len(gpr)-len(rxn_gpr)

#calculate the number of Promiscuous enzymes
#remove the reaction with isoenyzme
gem_rxn_nov1 = gem_rxn_nov[gem_rxn_nov['GPR'].str.contains('or', na=False)==False]
gem_rxn_nov2 = gem_rxn_nov1[gem_rxn_nov1['GPR'].str.contains('and', na=False)==False]
#remove the prot_ reaction
gem_rxn_nov3 = gem_rxn_nov2[gem_rxn_nov2['rxnID'].str.contains('prot_', na=False)==False]
promoiscuous = pd.Series(gem_rxn_nov3['GPR']).value_counts()
promoiscuous0 = promoiscuous[promoiscuous >=2]
len(promoiscuous0)-1