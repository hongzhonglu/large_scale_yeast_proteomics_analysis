# this module is mainly for ecModel simulation

from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *



# rerun the above step using a function
def DLecModelSimulate(model, dilution_rate):
   """
   This function is used to do simulation with ecModels using kcat value from deep learning.
   :param model: a ecModel
   :param dilution_rate: a dilution rate 0-0.42 /h

   :return: solution_f: fluxes datasets

   """

   dilutionrate = dilution_rate
   ecYeast = model
   if dilutionrate >= 0.4:
       ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27 * dilutionrate / 0.4, 0)  # this value is further rescaled by maximal growth rate at 0.42.
   else:
       ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)  # -230/0.55*0.4, this is rescaled by maximal growth rate.

   # refer to bioRxiv
   ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
              'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
              'oxygen exchange', 'EX_protein_pool']
   # find the related rxnID
   idx = []
   for name0 in ex_mets:
       print(name0)
       s = getRxnByReactionName(model=ecYeast, name=name0)
       if len(s) > 1:
           print("need check")
       elif len(s) == 1:
           idx.append(s[0])

   model_tmp = ecYeast.copy()
   model_tmp.reactions.get_by_id("r_1714").lower_bound = 0
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = -1000
   model_tmp.reactions.get_by_id(idx[0]).lower_bound = dilutionrate
   model_tmp.objective = {model_tmp.reactions.r_1714: 1}  # minimize the uptake of glucose
   solution2 = model_tmp.optimize()
   # then fix glucose uptake and minimize the protein pool
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = solution2.objective_value * 1.00001
   model_tmp.reactions.get_by_id(idx[9]).lower_bound = -1000
   model_tmp.objective = {model_tmp.reactions.EX_protein_pool: 1}  # minimize the usage of protein pools
   solution_f = model_tmp.optimize()
   solution_f.fluxes["EX_protein_pool"]
   solution_f.fluxes["r_1714"]
   return solution_f



# second ecYeast based om deep learning
# simulate the Crabtree effect
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
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


















# compare the predicted and measured protein abundances
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
solution3 = DLecModelSimulate(model=ecYeast, dilution_rate=0.42)
flux_max = solution3.fluxes
result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
result = result[result['rxnID'].str.contains("prot_")]
result['geneID'] = result['rxnID'].str.replace("prot_", "")


# input the proteomics under max growth rate
abundance_ex = pd.read_excel("data/proteomics/data_PNAS_2021.xlsx")
abundance_ex['g/gDW'] =(abundance_ex['replicate 1 (g gDW-1)']+ abundance_ex['replicate 2 (g gDW-1)']+ abundance_ex['replicate 3 (g gDW-1)'])/3
abundance_ex=abundance_ex[['Symbol','g/gDW']]
abundance_ex.columns = ['gene','g/gDW']
abundance_ex1 = splitAbundance(pro_df=abundance_ex)
# change the unit from g/gDW as mmol/gDW
# input the molecular weight
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
abundance_ex1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex1["gene"])
abundance_ex_check = abundance_ex1[abundance_ex1["MW_Kda"].isna()]
abundance_ex1=abundance_ex1[~abundance_ex1["MW_Kda"].isna()]
abundance_ex1["mmol/gDW"] = abundance_ex1["g/gDW"]/abundance_ex1["MW_Kda"]# #mmol/g biomass

result['pro_measured'] = singleMapping(abundance_ex1["mmol/gDW"],abundance_ex1["gene"],result['geneID'])
result = result[~result["pro_measured"].isna()]
result.to_excel("data/data_check.xlsx")


# change the protein abundance unit from mmol/gDW into protein copy/cell
coefficient1 = 7.8298e9
result_unify = result.copy()
result_unify["pro_measured"] = result['pro_measured']*coefficient1
result_unify["flux"] = result['flux']*coefficient1








# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

# method1 absolute protein abundance mmol protein/gDW
plt.figure()
sns.regplot(x=np.log10(result['pro_measured']), y=np.log10(result['flux']), fit_reg=False)
plt.xlim(-11, 0)
plt.ylim(-11, 0)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")


# method2 protein copy/cell
plt.figure()
sns.regplot(x=np.log10(result_unify['pro_measured']+1), y=np.log10(result_unify['flux']+1), fit_reg=False)
plt.xlim(-0.5, 7)
plt.ylim(-0.5, 7)
plt.xlabel("log10(Measured_protein_copy/cell + 1)")
plt.ylabel("log10(Predicted_protein_copy/cell +1)")

plt.figure()
sns.regplot(x=result_unify['pro_measured'], y=result_unify['flux'], fit_reg=False)
plt.xlabel("Measured_protein_copy/cell")
plt.ylabel("Predicted_protein_copy/cell")





# calculate the correlation coefficients - method1
# remove the proteins with zero
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)




# calculate the correlation coefficients - method2
corr, ss = pearsonr(np.log10(result_unify['pro_measured']+1), np.log10(result_unify['flux']+1))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)


result_unify1 = result_unify[result_unify['flux'] > 0]
result_unify1 = result_unify1[result_unify1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result_unify1['pro_measured']), np.log10(result_unify1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)
