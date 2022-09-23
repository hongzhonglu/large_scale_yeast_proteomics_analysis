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
dir2 = "data/ecYeast_DL_update_some_kcat.xml"
ecYeast = read_sbml_model(dir2)

# solve the model
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


result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)


# method2 protein copy/cell
plt.figure()
sns.regplot(x=np.log10(result_unify['pro_measured']+1), y=np.log10(result_unify['flux']+1), fit_reg=False)
plt.xlim(-0.5, 7)
plt.ylim(-0.5, 7)
plt.xlabel("log10(Measured_protein_copy/cell + 1)")
plt.ylabel("log10(Predicted_protein_copy/cell +1)")






# part 2
# classify metabolic genes based on organelles
# then check the correlation in each organelles
gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("prot_")]
gene_prot['geneID'] = gene_prot['rxnID'].str.replace("prot_", "")
ss = gene_prot[gene_prot['geneID'].str.contains("-")]
organelle_v = collectOrganelleTerm(type="volume")
organelle_m = collectOrganelleTerm(type="m")
gene_metabolic = gene_prot["geneID"].tolist()
compartment_in = organelle_v + organelle_m
m_gene_in_organelle = FingGenesForOrganelle(gene_set=gene_metabolic, compartment_list=compartment_in, compartment_type="organelle")


# loop
compartment_in0 = CompartmentInGEMs()

# Note: it shows that four genes from 'Golgi membrane' were not related to the core metabolic functions from the model, so the predicted protein abudance is zero.
result_unify = result_unify[result_unify['flux'] > 0]
correlation_all = []
for organelle_target in compartment_in0:
    print(organelle_target)
    gene_target = m_gene_in_organelle[organelle_target]
    result_unify_c = result_unify[result_unify["geneID"].isin(gene_target)]
    plt.figure()
    sns.regplot(x=np.log10(result_unify_c['pro_measured']), y=np.log10(result_unify_c['flux']), fit_reg=False)
    plt.xlim(-0.5, 7)
    plt.ylim(-0.5, 7)
    plt.xlabel("log10(Measured_protein_copy/cell)")
    plt.ylabel("log10(Predicted_protein_copy/cell)")
    plt.suptitle(organelle_target)
    # calculate the correlation coefficients - method2
    if result_unify_c.shape[0] >= 4:
        corr, ss = pearsonr(np.log10(result_unify_c['pro_measured']), np.log10(result_unify_c['flux']))
        print("Correlation coefficient:", corr)
        print("Correlation p_value:", ss)
        plt.savefig('result/figure/' + organelle_target + '_simulation_VS_measured.pdf', bbox_inches='tight')
    else:
        corr = None
    correlation_all.append(corr)

df = pd.DataFrame({"compartment":compartment_in0,"coefficient":correlation_all})
df = df.sort_values(by=['coefficient'], ascending=False)


# bar plot
plt.figure()
sns.set_style('darkgrid')
sns.barplot(x="compartment", y="coefficient", data=df, capsize=.2)
plt.xlabel("compartment", fontsize=12)
plt.ylabel("coefficient", fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()





# Note: It also shows that the predicted abundance of YJL167W is much higher than the measured one!
# 'endoplasmic reticulum' contains YJL167W
#getRxnByGene(ecYeast, "YJL167W")
#(sum(result_unify_c['flux'])-449824)/(sum(result_unify_c['pro_measured'])-84034)
#(sum(result_unify_c['flux']))/(sum(result_unify_c['pro_measured']))
# Note: endoplasmic reticulum membrane, predicted protein abundance for YNR016C and YGR060W is much higher than measured


