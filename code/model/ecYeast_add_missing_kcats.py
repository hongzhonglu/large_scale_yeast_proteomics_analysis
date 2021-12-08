# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import load_matlab_model
import os
import sys
import pprint
os.chdir('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code')
sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
pprint.pprint(sys.path)

# import self function
from src.mainFunction import *

dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"

ecYeast = load_matlab_model(dir1)

# check which protein has kcat and which protein has no kcat
gene_list = []
for x in ecYeast.genes:
    gene_list.append(x.id)

gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("draw_prot")]
gene_prot_list = gene_prot["GPR"].to_list()

gene_no_kinetic = list(set(gene_list)-set(gene_prot_list))

# TODO: there are 167 proteins with no kinetic information. Need additional check!
# how to get the missing kcat information