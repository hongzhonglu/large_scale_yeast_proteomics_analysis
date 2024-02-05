# Crabtree simulation

# import self function
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns

# load the model
ecYeast2 = getOrganelleConstraintGEM()
# check which organelle affect the crabtree?
# it shows that mitochondrion is main organelle to affect the output
organelle_pro_range = pd.read_excel("result/organelle_protein_abundance_range_rosemary.xlsx")
org1 = "mitochondrion"
constraint_name0 = org1 + '_constraint'
constraint_name0 = constraint_name0.replace(' ', '_')
print(constraint_name0)
compartment_info = organelle_pro_range[org1].tolist()


max_value = compartment_info[6] # adjust the upper bound
ecYeast2.constraints[constraint_name0].ub = max_value

# simulate crabtree effect!
# refer to bioRxiv
ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
           'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
           'oxygen exchange', 'EX_protein_pool']
# find the related rxnID
idx = []
for name0 in ex_mets:
    print(name0)
    s = getRxnByReactionName(model=ecYeast2, name=name0)
    if len(s) > 1:
        print("need check")
    elif len(s) == 1:
        idx.append(s[0])

# for the loop
dilutionrate = np.arange(0.05, 0.4, 0.05).tolist()
result = dict()
for k in range(len(dilutionrate)):
    print(k)
    model_tmp = ecYeast2
    solution3 = DLecModelSimulate(model=ecYeast2, dilution_rate=dilutionrate[k])
    fluxes = solution3.fluxes[idx]
    result[dilutionrate[k]] = list(fluxes)
# summarize the result
result_df = pd.DataFrame.from_dict(result)
result_df1 = result_df.transpose()
result_df1.columns = ex_mets
result_df1["D-glucose exchange"] = result_df1["D-glucose exchange"] * (-1)
result_df1["oxygen exchange"] = result_df1["oxygen exchange"] * (-1)
result_df2 = result_df1.drop('EX_protein_pool', 1)
# plot
plt.figure()
sns.lineplot(x='biomass pseudoreaction', y='value', hue='variable', style="variable",
             data=pd.melt(result_df2, ['biomass pseudoreaction']))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('result/figure/Cratree simulation based on ecModel_DLkcat.pdf', bbox_inches='tight')
