# This module is mainly used to build a pipeline to integrate structure information with models.


# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns


# second ecYeast based om deep learning
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("prot_")]
gene_prot['geneID'] = gene_prot['rxnID'].str.replace("prot_", "")




