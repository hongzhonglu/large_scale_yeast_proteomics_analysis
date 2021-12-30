# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1

from cobra.io import load_matlab_model

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"

# dir2 = "data/GEMs/ecYeastGEM_carl.mat" carl's model has no protein pool constraint for enzymes.

# it shows if the ecYeastGEM from ecModels was used, there will be some wierd things. It may due to some kcat of specific enzymes was not optimized!
# though this model has more gene numbers.
dir3 = "/Users/xluhon/Documents/GitHub/ecModels/ecYeastGEM/model/ecYeastGEM_batch.mat" # can't calculate the maximal growth at 0.424??

ecYeast = load_matlab_model(dir1)

# reaction annotation
# r_2111 growth 0.424
# r_1714_REV glucose uptake 18.79
# r_1992_REV oxygen uptake 2.92
# r_1672 co2 production 24.77
# r_1761 ethanol production 26.97
# r_1634 acetate secretion 0.566
# Glycerol secretion rate 1.469
# Pyruvate 0.163
# Succinate 0.0365
# growth 0.424


# protein_ratio 0.483


# simulation based on the minimization of protein abundances
# it should be noted that if the growth the over than 0.38, the model will not have right solution.
model = ecYeast.copy()
#gem_rxn_nov = produceRxnList(model)

model.reactions.get_by_id("r_1634").upper_bound = 0 # assume acetate is not produced!
model.reactions.get_by_id("r_2033").upper_bound = 0 # assume pyruvate is not produced!
model.reactions.get_by_id("r_1631").upper_bound = 0 # assume acetaldehyde is not produced!
model.reactions.get_by_id("r_1549").upper_bound = 0 # assume (R,R)-2,3-butanediol is not produced!
model.reactions.get_by_id("prot_pool_exchange").upper_bound = 0.10366*1.125  # assume (R,R)-2,3-butanediol is not produced!



# how to predict the max growth
model.objective = {model.reactions.r_2111: 1}
solution3 = model.optimize()

solution3.fluxes["r_1992_REV"]
solution3.fluxes["r_1672"]
solution3.fluxes["r_1761"]
solution3.fluxes["r_2111"]
solution3.fluxes["r_1714_REV"]
solution3.fluxes['prot_pool_exchange']



# next using the chemostat simulation
model.reactions.get_by_id("prot_pool_exchange").upper_bound = 0.10366*1.125  # assume (R,R)-2,3-butanediol is not produced!
i = 0.424
solution = chemostatSimulation(model0=model, D0=i)
solution.fluxes["r_1992_REV"]
solution.fluxes["r_1672"]
solution.fluxes["r_1761"]
solution.fluxes["r_2111"]
solution.fluxes["r_1714_REV"]
solution.fluxes['prot_pool_exchange']





# here compare the predicted proteomics and measured？
# plot the relation between the predict and measured proteins
flux_max = solution.fluxes
result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
result = result[result['rxnID'].str.contains("draw_prot")]
result['rxnID'] = result['rxnID'].str.replace("draw_prot_", "")
ID_map = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
result['geneID'] = singleMapping(ID_map['GeneName'], ID_map['Entry'], result['rxnID'])



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
result.to_excel("result/predicted_and_measured_proteomics.xlsx")

result=result[~result["pro_measured"].isna()]



# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

plt.figure()
sns.regplot(x=result['pro_measured'], y=result['flux'], fit_reg=False)
sns.regplot(x=np.log10(result['pro_measured']), y=np.log10(result['flux']), fit_reg=False)
plt.xlim(-10, -3)
plt.ylim(-10, -3)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")


# remove the proteins with zero
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)







## TO DO
# how to update the strain biomass in python?