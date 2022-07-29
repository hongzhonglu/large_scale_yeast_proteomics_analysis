# this module is mainly for ecModel simulation
# for the ecModels with organelle constraint, can it still predict Crabtree effect? In principle it should.



from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *


# second ecYeast based om deep learning
# simulate the Crabtree effect
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)


# using the manual curated ecYeast from deep learning
# in this version of model, we curate the kcat for some enzymes
ecYeast = read_sbml_model("data/ecYeast_DL_update_some_kcat.xml")




# refer to bioRxiv
ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
           'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange', 'oxygen exchange', 'EX_protein_pool']
# find the related rxnID
idx = []
for name0 in ex_mets:
    print(name0)
    s = getRxnByReactionName(model=ecYeast, name=name0)
    if len(s) > 1:
        print("need check")
    elif len(s) == 1:
        idx.append(s[0])

# for the loop
dilutionrate = np.arange(0.05, 0.42, 0.05).tolist()
result = dict()
for k in range(len(dilutionrate)):
    print(k)
    model_tmp = ecYeast.copy()
    solution3 = DLecModelSimulate(model=ecYeast, dilution_rate=dilutionrate[k])
    fluxes = solution3.fluxes[idx]
    result[dilutionrate[k]] = list(fluxes)


# summarize the result
result_df = pd.DataFrame.from_dict(result)
result_df1 = result_df.transpose()
result_df1.columns = ex_mets
result_df1["D-glucose exchange"] = result_df1["D-glucose exchange"]*(-1)
result_df1["oxygen exchange"] = result_df1["oxygen exchange"]*(-1)
result_df2 = result_df1.drop('EX_protein_pool', 1)
# plot
plt.figure()
sns.lineplot(x='biomass pseudoreaction', y='value', hue='variable', style="variable",
             data=pd.melt(result_df2, ['biomass pseudoreaction']))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('result/figure/Cratree simulation based on ecModel_DLkcat.pdf', bbox_inches='tight')
