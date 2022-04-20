# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import load_matlab_model
from cobra import Reaction, Metabolite

import sys
import pprint

pprint.pprint(sys.path)

# import self function
from src.mainFunction import *
from src.model_process import *
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
# check which reaction contains these genes with no kinetic parameters


GEM_subsystem = pd.read_csv('data/subsystem/Rxn_unique_subsystem_v2.tsv', sep="\t")
rxn_gene = getRXNgeneMapping(rxn0=GEM_subsystem['ID'], gpr0=GEM_subsystem['GENE ASSOCIATION'])
# analyze glucose
df = rxn_gene[rxn_gene['gene'].isin(gene_no_kinetic)]
GEM_select = GEM_subsystem[GEM_subsystem['ID'].isin(df['rxnID'])]
# it shows a lot of transport has no kcat data!



# try to add the glucose transporter
genes_select_glucose = getProteinForRxnGEM(rxnID=['r_1166'])
# input the molecular weight information
molecular_weight = pd.read_csv('data/sce_protein_weight.tsv', sep="\t")
molecular_weight_select = molecular_weight[molecular_weight['locus'].isin(genes_select_glucose)]












# a whole process to add the glucose transporter sectional area as constraints
model = ecYeast.copy()
# EC 2.5.1.32
MW1 = molecular_weight_select['proteins_molecular_weight'].mean()
kcat1 = 125 # s^(-1)


# firstly add the new metabolite
model.add_metabolites(Metabolite('prot_glc_trans',compartment='c', formula='', name=''))


# rewrite the reaction in transporting glucose through a pseudo glucose transporter
# as there are total about 20 glucose transporters, here we just use one pseudo transporter
dict1 = {model.metabolites.get_by_id('s_0565'): -1,
         model.metabolites.get_by_id('prot_glc_trans'): -1/(kcat1*3600),
         model.metabolites.get_by_id('s_0563'): 1,
        }


reaction = Reaction('pseudo_glc_trans')  # rxn ID
reaction.name = 'pseudo_glc_trans'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0  # This is the default
reaction.upper_bound = 1000  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict1)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])


# add the protein pool to the glucose transporter
dict2 = {model.metabolites.get_by_id('prot_pool'): -MW1/1000,
         model.metabolites.get_by_id('prot_glc_trans'):1
        }
reaction = Reaction('draw_pseudo_glc_trans')  # rxn ID
reaction.name = 'draw_pseudo_glc_trans'
reaction.subsystem = ''
reaction.lower_bound = 0  # This is the default
reaction.upper_bound = 1000  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict2)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])
# check the reaction
for r in model.reactions:
    print(r.id, r.name, r.reaction, r.compartments, sep="\t")





# try to simulate using the new models
model0 = ecYeastMinimalMedia(model)
model0.reactions.get_by_id("r_1634").upper_bound = 0 # assume acetate is not produced!
model0.reactions.get_by_id("r_2033").upper_bound = 0 # assume pyruvate is not produced!
model0.reactions.get_by_id("r_1631").upper_bound = 0 # assume acetaldehyde is not produced!
model0.reactions.get_by_id("r_1549").upper_bound = 0 # assume (R,R)-2,3-butanediol is not produced!
# turn off the original glucose transporter
model0.reactions.get_by_id("r_1166").bounds = (0, 0) # assume (R,R)-2,3-butanediol is not produced!




# here we can set the constraint for glucose transporter
coefficient1 = 6.5789e9
Sglucose_trans = 4.2 #4.2 # um^2 upper bound of sectional area occupied by glucose transporter
glucose_transporter_area = 20.344 # nm^2 average of glucose transporters sectional area in yeast
model0.reactions.get_by_id("draw_pseudo_glc_trans").upper_bound = 1e6*Sglucose_trans/(coefficient1*glucose_transporter_area)




# set growth
growth = 0.35
model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
# minimization glucose uptake rate
model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
model0.objective = {model0.reactions.r_1714_REV: -1}
solution2 = model0.optimize()
solution2.fluxes['pseudo_glc_trans']
solution2.fluxes['draw_pseudo_glc_trans']



GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
model0.objective = {model0.reactions.prot_pool_exchange: -1}
solution3 = model0.optimize()
solution3.fluxes['draw_pseudo_glc_trans']
solution3.fluxes['prot_pool_exchange']
solution3.fluxes["r_1761"]

