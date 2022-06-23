# This module is mainly used to build a pipeline to integrate structure information with models.


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



# all metabolic genes from ecGEMs
organelle_v = collectOrganelleTerm(type="volume")
organelle_m = collectOrganelleTerm(type="m")
gene_metabolic = gene_prot["geneID"].tolist()
compartment_in = organelle_v + organelle_m
m_gene_in_organelle = FingGenesForOrganelle(gene_set=gene_metabolic, compartment_list=compartment_in, compartment_type="organelle")


gene_plasma = m_gene_in_organelle['plasma membrane']
# find rxnid based on gene
rxn_select = gene_prot[gene_prot["geneID"].isin(gene_plasma)]["rxnID"].tolist()
rxn_select = [x.replace("-A", "_A") for x in rxn_select]
# generate the formula as the contraint
formula_list = ["ecYeast.reactions." + x + ".flux_expression" for x in rxn_select]
formula_one = " + ".join(formula_list)


# the following two seems consistent with each other
"""
same_flux = ecYeast.problem.Constraint(ecYeast.reactions.prot_YLR342W.flux_expression + ecYeast.reactions.prot_YPR159W.flux_expression + ecYeast.reactions.prot_YJL005W.flux_expression + ecYeast.reactions.prot_YCR024C_A.flux_expression + ecYeast.reactions.prot_YGL008C.flux_expression + ecYeast.reactions.prot_YPL036W.flux_expression + ecYeast.reactions.prot_YGR060W.flux_expression + ecYeast.reactions.prot_YNL192W.flux_expression + ecYeast.reactions.prot_YOR171C.flux_expression + ecYeast.reactions.prot_YOR317W.flux_expression + ecYeast.reactions.prot_YMR246W.flux_expression + ecYeast.reactions.prot_YJL100W.flux_expression + ecYeast.reactions.prot_YLR305C.flux_expression + ecYeast.reactions.prot_YDR208W.flux_expression + ecYeast.reactions.prot_YMR008C.flux_expression + ecYeast.reactions.prot_YLR020C.flux_expression + ecYeast.reactions.prot_YOR348C.flux_expression + ecYeast.reactions.prot_YCR010C.flux_expression + ecYeast.reactions.prot_YER056C.flux_expression + ecYeast.reactions.prot_YER060W.flux_expression + ecYeast.reactions.prot_YER060W_A.flux_expression + ecYeast.reactions.prot_YGL186C.flux_expression + ecYeast.reactions.prot_YJR152W.flux_expression + ecYeast.reactions.prot_YDR384C.flux_expression + ecYeast.reactions.prot_YGR121C.flux_expression + ecYeast.reactions.prot_YNL142W.flux_expression + ecYeast.reactions.prot_YGR065C.flux_expression + ecYeast.reactions.prot_YDR536W.flux_expression + ecYeast.reactions.prot_YLR081W.flux_expression + ecYeast.reactions.prot_YKL217W.flux_expression + ecYeast.reactions.prot_YDR342C.flux_expression + ecYeast.reactions.prot_YDR345C.flux_expression + ecYeast.reactions.prot_YHR092C.flux_expression + ecYeast.reactions.prot_YHR094C.flux_expression + ecYeast.reactions.prot_YMR011W.flux_expression + ecYeast.reactions.prot_YOR011W.flux_expression + ecYeast.reactions.prot_YCR098C.flux_expression + ecYeast.reactions.prot_YGL084C.flux_expression + ecYeast.reactions.prot_YLL043W.flux_expression + ecYeast.reactions.prot_YCL025C.flux_expression + ecYeast.reactions.prot_YKR039W.flux_expression + ecYeast.reactions.prot_YOL020W.flux_expression + ecYeast.reactions.prot_YDR497C.flux_expression + ecYeast.reactions.prot_YOL103W.flux_expression + ecYeast.reactions.prot_YMR319C.flux_expression + ecYeast.reactions.prot_YMR058W.flux_expression + ecYeast.reactions.prot_YBR068C.flux_expression + ecYeast.reactions.prot_YEL063C.flux_expression + ecYeast.reactions.prot_YBR069C.flux_expression + ecYeast.reactions.prot_YCR075C.flux_expression + ecYeast.reactions.prot_YGR191W.flux_expression + ecYeast.reactions.prot_YGR055W.flux_expression + ecYeast.reactions.prot_YCR028C.flux_expression + ecYeast.reactions.prot_YBR296C.flux_expression + ecYeast.reactions.prot_YCR037C.flux_expression + ecYeast.reactions.prot_YJL198W.flux_expression + ecYeast.reactions.prot_YJL129C.flux_expression + ecYeast.reactions.prot_YLL028W.flux_expression + ecYeast.reactions.prot_YOR273C.flux_expression + ecYeast.reactions.prot_YPL274W.flux_expression + ecYeast.reactions.prot_YLL061W.flux_expression + ecYeast.reactions.prot_YLR138W.flux_expression + ecYeast.reactions.prot_YHL016C.flux_expression + ecYeast.reactions.prot_YGR138C.flux_expression + ecYeast.reactions.prot_YPR156C.flux_expression + ecYeast.reactions.prot_YBR294W.flux_expression + ecYeast.reactions.prot_YLR092W.flux_expression + ecYeast.reactions.prot_YPL092W.flux_expression + ecYeast.reactions.prot_YLR237W.flux_expression + ecYeast.reactions.prot_YOR071C.flux_expression + ecYeast.reactions.prot_YBR021W.flux_expression + ecYeast.reactions.prot_YBL042C.flux_expression + ecYeast.reactions.prot_YLL052C.flux_expression + ecYeast.reactions.prot_YPR192W.flux_expression + ecYeast.reactions.prot_YNL065W.flux_expression + ecYeast.reactions.prot_YOR306C.flux_expression + ecYeast.reactions.prot_YMR162C.flux_expression + ecYeast.reactions.prot_YBR295W.flux_expression + ecYeast.reactions.prot_YML125C.flux_expression + ecYeast.reactions.prot_YDR038C.flux_expression + ecYeast.reactions.prot_YDR039C.flux_expression + ecYeast.reactions.prot_YDR040C.flux_expression + ecYeast.reactions.prot_YJR040W.flux_expression + ecYeast.reactions.prot_YKL220C.flux_expression + ecYeast.reactions.prot_YLR214W.flux_expression + ecYeast.reactions.prot_YNR060W.flux_expression + ecYeast.reactions.prot_YOL152W.flux_expression + ecYeast.reactions.prot_YOR381W.flux_expression + ecYeast.reactions.prot_YKR093W.flux_expression + ecYeast.reactions.prot_YDR093W.flux_expression + ecYeast.reactions.prot_YOL122C.flux_expression + ecYeast.reactions.prot_YNL275W.flux_expression + ecYeast.reactions.prot_YLR130C.flux_expression + ecYeast.reactions.prot_YOL130W.flux_expression,
    lb=0.00002,
    ub=0.00006)
ecYeast.add_cons_vars(same_flux)"""

same_flux = ecYeast.problem.Constraint(
    eval(formula_one),
    lb=0.00002,
    ub=0.00006)
#ecYeast.add_cons_vars(same_flux)
ecYeast.constraints.append(same_flux) # this is a right way to add additional constraint


# try to simulation with this new constraint
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

















