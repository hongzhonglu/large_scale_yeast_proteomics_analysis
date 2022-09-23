# this module is mainly for ecModel simulation

from cobra.io import read_sbml_model
import cobra

# import self function
from src.mainFunction import *
from src.protein_process import *
from src.model_process import *

# compare the predicted and measured protein abundances
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)

# update the model
# manual curation 1
# rxnID = getRxnByGene(model=ecYeast, gene0='YNR016C')
target_gene0 = 'YNR016C'
rxnID0 = 'r_0109'
kcat_m0 = 8.5*3600 # unit is /h
ecYeast = updateEcGEMkcat(ecGEM=ecYeast, target_gene=target_gene0, rxnID=rxnID0, kcat_m=kcat_m0)

# manual curation 2
rxnID = getRxnByGene(model=ecYeast, gene0='YGR060W')
# using a loop
for rxn0 in rxnID:
    target_gene0 = 'YGR060W'
    kcat_m0 = 1/6.11002831284031e-05
    ecYeast = updateEcGEMkcat(ecGEM=ecYeast, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)

# save the model
cobra.io.write_sbml_model(ecYeast, "data/ecYeast_DL_update_some_kcat.xml")








