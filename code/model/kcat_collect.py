import pandas as pd
from src.mainFunction import *
kcat_DP = pd.read_table('data/kcat_from_deep_learning/Saccharomyces_cerevisiae_PredictionResults.txt')


gene_rxn = splitAndCombine(gene=kcat_DP['genes'], rxn=kcat_DP['# rxnID'], sep0=";")
# check the unique genes
all_gene = list(set(gene_rxn['V2'].to_list()))


















