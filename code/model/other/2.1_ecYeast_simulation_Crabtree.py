# this module is mainly for crabtree simulation using ecModel from GECKO
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
from src.model_process import *
from src.mainFunction import *

def Crabtree_Simulate(ecModel_batch, growth_list):
    model = ecModel_batch.copy()
    growth = growth_list
    glucose_uptake = []
    oxgen_uptake = []
    co2_production = []
    ethanol_production = []
    acetate_production = []
    protein_pool = []
    model.reactions.get_by_id("r_1634").upper_bound = 0  # assume acetate is not produced!
    model.reactions.get_by_id("r_2033").upper_bound = 0  # assume pyruvate is not produced!
    model.reactions.get_by_id("r_1631").upper_bound = 0  # assume acetaldehyde is not produced!
    model.reactions.get_by_id("r_1549").upper_bound = 0  # assume (R,R)-2,3-butanediol is not produced!
    gem_rxn_nov = produceRxnList(model)

    for i in growth:
    #    if i < 0.43:
           # print(i)
        string0 = "D=" + str(i)
        solution = chemostatSimulation(model0=model, D0=i)
        oxgen_uptake.append(solution.fluxes["r_1992_REV"])
        co2_production.append(solution.fluxes["r_1672"])
        ethanol_production.append(solution.fluxes["r_1761"])
        acetate_production.append(solution.fluxes["r_1634"])
        glucose_uptake.append(solution.fluxes["r_1714_REV"])
        protein_pool.append(solution.fluxes['prot_pool_exchange'])
        gem_rxn_nov[string0] = list(solution.fluxes)

    #    else:
    #        break

    # plot the relation between the glucose uptake rate and growth rate
    # figure 1
    plt.figure()
    plt.plot(growth, glucose_uptake, marker='.', label='Glucose')
    plt.plot(growth, oxgen_uptake, marker='.', label='O2')
    plt.plot(growth, co2_production, marker='.', label='CO2')
    plt.plot(growth, ethanol_production, marker='.', label='ethanol')
    plt.plot(growth, acetate_production, marker='.', label='acetate')
    plt.xlabel('Growth rate (/h)')
    plt.ylabel('rate (mmol/gDW.h)')
    plt.legend(loc='upper left')
    plt.ylim(0, 45)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("result/figure/crabtree_effect_simulation.pdf", bbox_inches='tight')

    # from the above simulation it could find the secretion of acetate
    #gem_rxn_nov.to_excel("result/flux_check.xlsx")
    # figure 2
    plt.figure()
    plt.plot(growth, protein_pool, marker='.', label='protein_pool')
    plt.xlabel('Growth rate (/h)')
    plt.ylabel('rate (mmol/gDW)')
    plt.legend(loc='upper left')
    plt.show()
    return gem_rxn_nov

dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"
#dir3 = "/Users/xluhon/Documents/GitHub/ecModels/ecYeastGEM/model/ecYeastGEM_batch.mat" # can't calculate the maximal growth at 0.424??
ecYeast = load_matlab_model(dir1)
# reaction annotation
# r_2111 growth
# r_1714_REV glucose uptake
# r_1992_REV oxygen uptake
# r_1672 co2 production
# r_1761 ethanol production
# r_1634 acetate secretion

# simulation based on the minimization of protein abundances
# it should be noted that if the growth the over than 0.38, the model will not have right solution.
model = ecYeast.copy()
growth0 = [0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.379]
fluxes = Crabtree_Simulate(ecModel_batch=model, growth_list=growth0)

# crabtree simulation for Jianye datasets
growth2 = [0.027, 0.044, 0.102, 0.152, 0.214, 0.254, 0.284, 0.334, 0.379]
fluxes = Crabtree_Simulate(ecModel_batch=model, growth_list=growth2)























# other scripts
# combine the predicted proteomics and measured
all_dilution = ["D=" + str(i) for i in growth2]
omcis_jianye = pd.read_csv("data/proteomics/Omics_from_Jianye.csv")
columns0 = list(omcis_jianye.columns)
columns0 = [x for x in columns0 if "RNA" not in x]
columns0 = [x for x in columns0 if "XIA" not in x]
omcis_jianye0 = omcis_jianye[columns0]
new_columns0 = ['Accession','Gene'] + [x + "_M" for x in all_dilution]
omcis_jianye0.columns = new_columns0
omcis_jianye1 = omcis_jianye0[['Accession','Gene']]
for x in new_columns0:
    if "_M" in x:
        print(x)
        ss1 = omcis_jianye0[x]*1e-09
        omcis_jianye1[x] = list(ss1)
result = fluxes[fluxes['rxnID'].str.contains("draw_prot")]
result['rxnID'] = result['rxnID'].str.replace("draw_prot_", "")
result0 = result[['GPR', 'rxnID'] + all_dilution]
# combine the predicted and measured omics
df_combine = pd.merge(left=result0, right=omcis_jianye1, left_on=['rxnID'], right_on = ['Accession'], how="left")
df_combine.to_excel("result/data_combine_xia_dataset.xlsx")


# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr


# check all
for x in all_dilution:
    p1 = x
    m1 = p1 + '_M'
    new_df = df_combine[[p1, m1]]
    new_df = new_df[~new_df[m1].isna()]
    plt.figure()
    sns.regplot(x=np.log10(new_df[m1]), y=np.log10(new_df[p1]), fit_reg=False)
    plt.xlim(-10, -2)
    plt.ylim(-10, -2)
    plt.xlabel("log10(Measured_protein_level)")
    plt.ylabel("log10(Predicted_protein_usage)")
    new_df1 = new_df[new_df[m1] > 0]
    new_df1 = new_df1[new_df1[p1] > 0]
    corr, ss = pearsonr(np.log10(new_df1[m1]), np.log10(new_df1[p1]))
    print("Correlation coefficient:", corr)
    print("Correlation p_value:", ss)
    plt.title(x + " Coefficient: " + str(corr))


