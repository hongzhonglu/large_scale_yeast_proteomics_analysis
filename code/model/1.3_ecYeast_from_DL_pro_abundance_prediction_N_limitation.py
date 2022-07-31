# this module is mainly for ecModel simulation

from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import sys
import matplotlib.pyplot as plt
import seaborn as sns
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




# solve the model
solution3 = DLecModelSimulate(model=ecYeast, dilution_rate=0.35)

flux_max = solution3.fluxes
result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
result = result[result['rxnID'].str.contains("prot_")]
result['geneID'] = result['rxnID'].str.replace("prot_", "")




# input the measured protein abundances
omics_tao2 = pd.read_excel("data/proteomics/Omics_from_tao_scale.xlsx")
condition = ['gene','prot.19', 'prot.20', 'prot.21']
omics_select = omics_tao2[condition]
omics_select['average'] = omics_select.drop('gene', axis=1).apply(lambda x: x.mean(), axis=1)

# compare the predict with measured
result['pro_measured'] = singleMapping(omics_select["average"], omics_select["gene"], result['geneID'])
result = result[~result["pro_measured"].isna()]

# method1 absolute protein abundance mmol protein/gDW
plt.figure()
sns.regplot(x=np.log10(result['pro_measured']), y=np.log10(result['flux']), fit_reg=False)
plt.xlim(-11, 0)
plt.ylim(-11, 0)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")

from scipy.stats import pearsonr
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)



# method2 calculate the RMSE
from sklearn.metrics import mean_squared_error
import math

y_actual = [1, 2, 3, 4, 5]
y_predicted = [1.6, 2.5, 2.9, 3, 4.1]
MSE = mean_squared_error(y_actual, y_predicted)
RMSE = math.sqrt(MSE)
print("Root Mean Square Error:\n")
print(RMSE)