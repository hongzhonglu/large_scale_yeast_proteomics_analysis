# Compare the predicted and measured protein abundance under N limitaiton

# import self function
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns

# load the model
ecYeast2 = getOrganelleConstraintGEM()

# just initial compare the predicted protein abundance and the total abundances
solution3 = DLecModelSimulate(model=ecYeast2, dilution_rate=0.35)
flux_max = solution3.fluxes
result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
result = result[result['rxnID'].str.contains("prot_")]
result['geneID'] = result['rxnID'].str.replace("prot_", "")

# input the measured protein abundances
omics_tao2 = pd.read_excel("data/proteomics/Omics_from_tao_scale.xlsx")
condition = ['gene','prot.19', 'prot.20', 'prot.21'] # growth rate 0.35/h
#condition = ['gene','prot.7', 'prot.8', 'prot.9'] # growth rate 0.10/h
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



# check the correlation and RMSE
# Correlation
from scipy.stats import pearsonr
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)
# RMSE
from sklearn.metrics import mean_squared_error
import math
y_actual = np.log10(result1['pro_measured'])
y_predicted = np.log10(result1['flux'])
MSE = mean_squared_error(y_actual, y_predicted)
RMSE = math.sqrt(MSE)
print("Root Mean Square Error:\n")
print(RMSE)
