# This pipeline could be used to evaluate the absolute protein abundance from each organelle affect the growth rate

# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns


# second ecYeast based om deep learning
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("prot_")]
gene_prot['geneID'] = gene_prot['rxnID'].str.replace("prot_", "")
ss = gene_prot[gene_prot['geneID'].str.contains("-")]

for rxn in ecYeast.reactions:
    if "-A" in rxn.id:
        print(rxn.id)
        ecYeast.reactions.get_by_id(rxn.id).id = rxn.id.replace("-A", "_A")


# using the manual curated ecYeast from deep learning
# in this version of model, we curate the kcat for some enzymes
ecYeast = read_sbml_model("data/ecYeast_DL_update_some_kcat.xml")

for rxn in ecYeast.reactions:
    if "-A" in rxn.id:
        print(rxn.id)
        ecYeast.reactions.get_by_id(rxn.id).id = rxn.id.replace("-A", "_A")



# generate the general formula as the constraint
# all metabolic genes from ecGEMs
organelle_v = collectOrganelleTerm(type="volume")
organelle_m = collectOrganelleTerm(type="m")
gene_metabolic = gene_prot["geneID"].tolist()
compartment_in = organelle_v + organelle_m
m_gene_in_organelle = FingGenesForOrganelle(gene_set=gene_metabolic, compartment_list=compartment_in, compartment_type="organelle")



# The following script is mainly used to check the organelle total abundance affect the cellular growth rate
# it is found some organlle protein total abudance have no affect while others have?
def AddOrgConstraint(ecModel, flux_expression):
    model_tmp = ecModel.copy()
    lower = 0
    upper = 1000
    model_tmp.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)
    same_flux = model_tmp.problem.Constraint(eval(flux_expression), lb=lower, ub=upper, name=constraint_name)
    model_tmp.add_cons_vars(same_flux)
    return model_tmp



# Add constraints in loop using the following scripts
# input the dataset information
organelle_pro_range = pd.read_excel("result/organelle_protein_abundance_range_rosemary.xlsx")
organelle_v0 = ['mitochondrion', 'nucleus', 'cytosol',
 'endoplasmic reticulum', 'lipid droplet', 'fungal-type vacuole',
 'peroxisome', 'Golgi apparatus']
organelle_m0 = ['fungal-type vacuole membrane',
 'plasma membrane',
 'mitochondrial outer membrane',
 'endoplasmic reticulum membrane',
 'mitochondrial inner membrane',
 'Golgi membrane',
 'peroxisomal membrane', #only with one metabolic gene from ecYeast
 'nuclear membrane'] #only with one metabolic gene from ecYeast
compartment_in0 = organelle_v0 + organelle_m0


for org in compartment_in0:
    print(org)
    organelle_target = org
    gene_target = m_gene_in_organelle[organelle_target]
    rxn_select = gene_prot[gene_prot["geneID"].isin(gene_target)]["rxnID"].tolist()
    rxn_select = [x.replace("-A", "_A") for x in rxn_select]
    formula_list = ["model_tmp.reactions." + x + ".flux_expression" for x in rxn_select]
    formula_one = " + ".join(formula_list)
    constraint_name = org + '_constraint'
    constraint_name = constraint_name.replace(' ','_')
    ecYeast = AddOrgConstraint(ecModel=ecYeast, flux_expression=formula_one)

# check the growth
objective = ecYeast.problem.Objective(ecYeast.reactions.r_4041.flux_expression, direction='max') # biomass
ecYeast.objective = objective
solution2 = ecYeast.optimize()
print("Max growth:", solution2.objective_value)





# reset the constraints????
ecYeast2 = ecYeast.copy()
for org in compartment_in0:
    constraint_name = org + '_constraint'
    constraint_name = constraint_name.replace(' ','_')
    print(constraint_name)
    compartment_info = organelle_pro_range[org].tolist()
    ecYeast2.constraints[constraint_name].ub = 1000
    ecYeast2.constraints[constraint_name].lb = 0
    if org == 'endoplasmic reticulum':
        max_value = compartment_info[6] #* 0.44 * 2.5
        ecYeast2.constraints[constraint_name].ub = max_value
        min_value = compartment_info[3]  # minimum  value
    elif org=='endoplasmic reticulum membrane':
        max_value = compartment_info[7] * 2.5  # 9 for origninal ecYeast2 from DL

        ecYeast2.constraints[constraint_name].ub = max_value
    else:
        min_value = compartment_info[3]  # minimum  value
        max_value = compartment_info[6]  # max value
        ecYeast2.constraints[constraint_name].lb = min_value
        ecYeast2.constraints[constraint_name].ub = max_value


# check the max growth
objective = ecYeast2.problem.Objective(ecYeast2.reactions.r_4041.flux_expression, direction='max') # biomass
ecYeast2.objective = objective
solution2 = ecYeast2.optimize()
print("Max growth:", solution2.objective_value)
## based on the above model with organelle constraint,
## next we use the general function to solve the model
# solve the model
solution3 = DLecModelSimulate(model=ecYeast2, dilution_rate=0.35)
#solution3 = DLecModelSimulate(model=ecYeast2, dilution_rate=0.1)




# just initial compare the predicted protein abundance and the total abundances
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

