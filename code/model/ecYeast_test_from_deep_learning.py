# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite

import sys
import pprint

pprint.pprint(sys.path)

# import self function
from src.mainFunction import *
from src.model_process import *

# some test functions
# Note: these following two functions are just for test analysis
def ecYeastMinimalMedia_No_reverse(model):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """
    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if '_REV' in x and 'exchange' in y]
    # first block any uptake
    for i, x in enumerate(exchange_rxn):
        rxn0 = exchange_rxn[i]
        #print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    #Allow uptake of essential components
    model.reactions.get_by_id("r_1654").lower_bound = -10000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861").lower_bound = -10000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100").lower_bound = -10000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992").lower_bound = -10000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005").lower_bound = -10000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060").lower_bound = -10000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832").lower_bound = -10000 #H+ exchange (reversible)
    return model

def chemostatSimulation_No_reverse(model0, D0):
    """
    This funcion is used to simulate the chemostat growth of yeast
    Actually this function is general to solve the ecGEMs
    :param model0: a ecGEMs
    :param D0: a growth rate
    :return: solution of fluxes
    """
    growth = D0
    with model0:
        model0 = ecYeastMinimalMedia_No_reverse(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714").bounds = (-1000, 100)  # open the glucose
        model0.objective = {model0.reactions.r_1714: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714").bounds = (GR * 1.001, GR)
        model0.objective = {model0.reactions.EX_protein_pool: -1}
        solution3 = model0.optimize()
    return solution3


# second ecYeast based om deep learning
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
# some initial constraints adjustment
ecYeast.reactions.get_by_id("r_4527").bounds = (0,0)  # ammonium exchange (reversible)
ecYeast.reactions.get_by_id("r_4538").bounds = (0,0)  # ammonium exchange (reversible)
ecYeast.reactions.get_by_id("r_4502").bounds = (0,0)  # ammonium exchange (reversible)
ecYeast.reactions.get_by_id("r_4504").bounds = (0,0)  # ammonium exchange (reversible)
ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)  # -230/0.55*0.4

# refer to bioRxiv
ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
           'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
           'oxygen exchange', 'EX_protein_pool']

# find the related rxnID
def getRxnByReactionName(model, name):
    s = []
    for rxn in model.reactions:
        if name == rxn.name:
            #print(rxn.id)
            s.append(rxn.id)
    return s


idx = []
for name0 in ex_mets:
    print(name0)
    s = getRxnByReactionName(model=ecYeast, name=name0)
    if len(s) > 1:
        print("need check")
    elif len(s) == 1:
        idx.append(s[0])


dilutionrate = np.arange(0.05, 0.42, 0.05).tolist()

result = dict()
for k in range(len(dilutionrate)):
    print(k)
    model_tmp = ecYeast.copy()
    model_tmp.reactions.get_by_id("r_1714").lower_bound = 0
    model_tmp.reactions.get_by_id(idx[1]).lower_bound = -1000
    model_tmp.reactions.get_by_id(idx[0]).lower_bound = dilutionrate[k]
    model_tmp.objective = {model_tmp.reactions.r_1714: 1}
    solution2 = model_tmp.optimize()
    solution2.fluxes["r_1714"]
    # then fix glucose uptake and minimize the protein pool
    model_tmp.reactions.get_by_id(idx[1]).lower_bound = solution2.objective_value*1.00001
    model_tmp.reactions.get_by_id(idx[9]).lower_bound = -1000
    model_tmp.objective = {model_tmp.reactions.EX_protein_pool: 1}
    solution3 = model_tmp.optimize()
    solution3.fluxes["EX_protein_pool"]
    solution3.fluxes["r_1714"]
    fluxes = solution3.fluxes[idx]
    result[dilutionrate[k]] = list(fluxes)

result_df = pd.DataFrame.from_dict(result)

result_df1 = result_df.transpose()
result_df1.columns = ex_mets



# plot
import matplotlib.pyplot as plt
import seaborn as sns
plt.figure()
sns.lineplot(x='biomass pseudoreaction', y='value', hue='variable', style="variable",
             data=pd.melt(result_df1, ['biomass pseudoreaction']))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)













# self script
i = 0.35
solution = chemostatSimulation_No_reverse(model0=ecYeast, D0=i)


# check which protein has kcat and which protein has no kcat
gene_list = []
for x in ecYeast.genes:
    gene_list.append(x.id)

gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("prot_")]
gene_prot_list = gene_prot["rxnID"].str.replace("prot_", "")
gene_no_kinetic = list(set(gene_list)-set(gene_prot_list))


# it means that all gene have kinetic information, so how to use the kineitc information and the models??
gem_rxn_nov["fluxes"] = list(solution2.fluxes)
gem_rxn_nov_exchange = gem_rxn_nov[gem_rxn_nov['name'].str.contains("exchange")]










