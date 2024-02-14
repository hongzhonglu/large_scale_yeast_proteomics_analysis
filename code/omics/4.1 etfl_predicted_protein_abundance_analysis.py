# Note: in this module, we will analyze the predicted protein abundance from etfl
# import self function
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import numpy

# part1
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx") # the unit the g/gDW???? should be wrong

# input the predicted abundance and combine it together
p1 = pd.read_excel("data/etfl_output/cefl_chemostat_enzymeUsage.xlsx")
p2 = pd.read_excel("data/etfl_output/cefl_batch_enzymeUsage.xlsx")
p1.rename(columns={p1.columns[0]: "all_gene"}, inplace = True)
p2.rename(columns={p2.columns[0]: "all_gene"}, inplace = True)
predicted = pd.merge(p1, p2, left_on="all_gene", right_on="all_gene", how="left")
predicted = predicted[predicted["all_gene"] != "dummy_peptide"]
predicted.to_excel("data/etfl_output/etfl_predicted_enzyme_usage.xlsx")

# condition_list
all_condition = list(predicted.columns)
measured = omics_combine_auto[all_condition]
# select the values
measured_simple0 = predicted[["all_gene"]]
measured_simple0 = pd.merge(measured_simple0, measured, left_on="all_gene", right_on="all_gene", how="left")

# analyze used and unused proteins
all_condition0 = all_condition[1:]
all_gene_predicted_in_use = []
number_of_used_gene = []
all_gene_list = []

for xx in all_condition0:
    print(xx)
    x0 = measured_simple0[xx]
    y0 = predicted[xx]
    combine_df = pd.DataFrame({"all_gene": predicted["all_gene"], "measured":x0, "predicted":y0})
    combine_df_used = combine_df[combine_df["predicted"] > 0]
    all_gene_predicted_in_use = all_gene_predicted_in_use + combine_df_used["all_gene"].tolist()
    number_of_used_gene.append(len(combine_df_used["all_gene"].tolist()))
    all_gene_list.append(combine_df_used["all_gene"].tolist())




# it is found that total 609 enzymes used under different nitrogen and carbon sources
# In average, 450 enzymes used under different conditions
# core enzymes: 196 core enzymes
all_gene_predicted_in_use = list(set(all_gene_predicted_in_use))
numpy.mean(number_of_used_gene)
core_enzyme = list(set.intersection(*[set(list_) for list_ in all_gene_list]))
",".join(core_enzyme)




# firstly calculate the correlation as a whole?
# plot
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr



all_corr = []
for xx in all_condition0:
    print(xx)
    x0 = measured_simple0[xx]
    y0 = predicted[xx]
    # remove the proteins with zero
    combine_df = pd.DataFrame({"all_gene": predicted["all_gene"], "measured": x0, "predicted": y0})
    result1 = combine_df[combine_df['measured'] > 0]  # for the log calculation, the value must be larger than zero
    result1 = result1[result1['predicted'] > 0]
    corr, ss = pearsonr(np.log10(result1['measured']), np.log10(result1['predicted']))
    print("Correlation coefficient:", corr)
    print("Correlation p_value:", ss)
    all_corr.append(corr)

corr_df = pd.DataFrame({"condition": all_condition0, "correlation": all_corr})
corr_df.to_excel("result/all_correlation_between_measured_protein_abundance_and_predicted.xlsx")




x0 = measured_simple0["D=0.027_M"]
y0 = predicted["D=0.027_M"]
# remove the proteins with zero
combine_df = pd.DataFrame({"all_gene":predicted["all_gene"], "measured":x0, "predicted":y0})
result1 = combine_df[combine_df['measured'] > 0] # for the log calculation, the value must be larger than zero
result1 = result1[result1['predicted'] > 0]
corr, ss = pearsonr(np.log10(result1['measured']), np.log10(result1['predicted']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)
plt.figure()
sns.regplot(x=np.log10(result1["measured"]), y=np.log10(result1["predicted"]), fit_reg=False)
plt.xlim(-10, -3)
plt.ylim(-10, -3)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")



