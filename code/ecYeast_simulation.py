# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import read_sbml_model, load_matlab_model
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pprint
os.chdir('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code')
sys.path.append(r"/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/code")
pprint.pprint(sys.path)

# import self function
from mainFunction import *

# function copy from strain_design repo
def ecYeastMinimalMedia(model):
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
    model.reactions.get_by_id("r_1654_REV").upper_bound = 10000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861_REV").upper_bound = 10000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100_REV").upper_bound = 10000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992_REV").upper_bound = 10000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005_REV").upper_bound = 10000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060_REV").upper_bound = 10000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832_REV").upper_bound = 10000 #H+ exchange (reversible)
    return model

def frange(start, stop, step):
    """
    This function is like range, step can be float value, like 0.1
    :param start:
    :param stop:
    :param step:
    :return:
    """
    i = start
    while i < stop:
        yield i
        i += step

# chemostat simulations
def chemostatSimulation(model0, D0):
    growth = D0
    with model0:
        model0 = ecYeastMinimalMedia(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    return solution3
dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"

ecYeast = load_matlab_model(dir1)

# reaction annotation
# r_2111 growth
# r_1714_REV glucose uptake
# r_1992_REV oxygen uptake
# r_1672 co2 production
# r_1761 ethanol production
# r_1634 acetate secretion


# simulation based on the minimization of protein abundances
model = ecYeast.copy()
growth = [0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36]
glucose_uptake = []
oxgen_uptake = []
co2_production = []
ethanol_production = []
acetate_production = []
protein_pool = []
model.reactions.get_by_id("r_1634").upper_bound = 0 # assume acetate is not produced!
model.reactions.get_by_id("r_2033").upper_bound = 0 # assume pyruvate is not produced!
model.reactions.get_by_id("r_1631").upper_bound = 0 # assume acetaldehyde is not produced!
model.reactions.get_by_id("r_1549").upper_bound = 0 # assume (R,R)-2,3-butanediol is not produced!

gem_rxn_nov = produceRxnList(model)

for i in growth:
    print(i)
    if i < 0.36:
        string0 = "D=" + str(i)
        solution = chemostatSimulation(model0=model, D0=i)
        oxgen_uptake.append(solution.fluxes["r_1992_REV"])
        co2_production.append(solution.fluxes["r_1672"])
        ethanol_production.append(solution.fluxes["r_1761"])
        acetate_production.append(solution.fluxes["r_1634"])
        growth.append(solution.fluxes["r_2111"])
        glucose_uptake.append(solution.fluxes["r_1714_REV"])
        protein_pool.append(solution.fluxes['prot_pool_exchange'])
        gem_rxn_nov[string0]=list(solution.fluxes)
    else:
        break

# plot the relation between the glucose uptake rate and growth rate
growth = growth[0:13]
plt.figure()
plt.plot(growth, glucose_uptake, marker='.', label='Glucose')
plt.plot(growth, oxgen_uptake, marker='.',label='O2')
plt.plot(growth, co2_production, marker='.',label='CO2')
plt.plot(growth, ethanol_production, marker='.',label='ethanol')
plt.plot(growth, acetate_production, marker='.',label='acetate')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.ylim(0, 25)
plt.show()

# from the above simulation it could find the secretion of acetate
gem_rxn_nov.to_excel("../result/flux_check.xlsx")

plt.figure()
plt.plot(growth, protein_pool, marker='.',label='protein_pool')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.show()
# from the above simulation it could find the secretion of acetate











# FBA analysis
model = ecYeast.copy()
model.objective = 'r_4041' # r_2111
model.reactions.get_by_id("r_1714_REV").bounds = (0, 5)
solution1 = model.optimize()
model.summary()

#combine the fluxes with rxn information
flux0 = pd.DataFrame(solution1.fluxes)
flux1 = flux0.rename_axis("rxnID").reset_index()
gem_rxn_nov = produceRxnList(model)
gem_rxn_nov["flux"] = flux1["fluxes"]
gem_rxn_nov["rxnID2"] = flux1["rxnID"]
gem_rxn_nov.to_excel("../result/flux_manual_check.xlsx")

# simulate the growth
"""analysis based on general FBA"""
# Open the glucose uptake
glucose_uptake = list(frange(0, 15, 0.2))
growth = []
oxgen_uptake = []
co2_production = []
ethanol_production = []
acetate_production = []
protein_pool = []
model0 = ecYeast.copy()
model0 = ecYeastMinimalMedia(model0)
model0.objective = 'r_2111'
model0.reactions.get_by_id("r_2111").lower_bound = 0
model0.reactions.get_by_id("r_1634").upper_bound = 0 # assume acetate is not produced!

for i in glucose_uptake:
    model0.reactions.get_by_id("r_1714_REV").bounds = (0, i) #D-glucose exchange (reversible)
    solution = model0.optimize()
    oxgen_uptake.append(solution.fluxes["r_1992_REV"])
    co2_production.append(solution.fluxes["r_1672"])
    ethanol_production.append(solution.fluxes["r_1761"])
    acetate_production.append(solution.fluxes["r_1634"])
    growth.append(solution.fluxes["r_2111"])
    protein_pool.append(solution.fluxes['prot_pool_exchange'])

# plot the relation between the glucose uptake rate and growth rate
plt.figure()
plt.plot(growth, glucose_uptake, marker='.', label='Glucose')
plt.plot(growth, oxgen_uptake, marker='.',label='O2')
plt.plot(growth, co2_production, marker='.',label='CO2')
plt.plot(growth, ethanol_production, marker='.',label='ethanol')
plt.plot(growth, acetate_production, marker='.',label='acetate')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.show()
# from the above simulation it could find the secretion of acetate




# analyze another models based on deep learning
# can't be input through cobrapy
dir2 = "../data/deep_learning_data/emodel_Saccharomyces_cerevisiae_Posterior_mean.mat"
ecYeast_DL = load_matlab_model(dir2)