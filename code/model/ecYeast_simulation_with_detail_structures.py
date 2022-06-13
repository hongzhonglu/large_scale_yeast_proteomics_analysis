# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import load_matlab_model
import matplotlib.pyplot as plt
import os
import sys
import pprint
#os.chdir('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code')
#sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
#pprint.pprint(sys.path)

# import self function
from src.mainFunction import *
from src.model_process import *


dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"
ecYeast = load_matlab_model(dir1)

# reaction annotation
# r_2111 growth
# r_1714_REV glucose uptake
# r_1992_REV oxygen uptake
# r_1672 co2 production
# r_1761 ethanol production
# r_1634 acetate secretion



#gem_rxn_nov = produceRxnList(ecYeast)
#gem_rxn_nov.to_excel('data/gem_rxn_nov.xlsx')

gem_rxn_nov = pd.read_excel('data/gem_rxn_nov.xlsx')
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("draw_prot")]
test_rxn = gene_prot["name"].to_list()

#then get the protein volume information
#input the protein volume datasets
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
gene_prot["Volume"] = singleMapping(pro_size['Total_Volume'], pro_size['locus'], gene_prot['GPR'])
gene_prot["section_area"] = singleMapping(pro_size['section_area_new'], pro_size['locus'], gene_prot['GPR'])


# first analyze the proteins for specific rxn
# now the model has no kinetic information for the glucose
genes_select_glucose = getProteinForRxnGEM(rxnID=['r_1166'])
gene_prot0 = gene_prot[gene_prot['GPR'].isin(genes_select_glucose)]
