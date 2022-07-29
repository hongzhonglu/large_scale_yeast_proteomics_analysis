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

from scipy.stats import pearsonr
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

plt.figure()
sns.regplot(x=result_unify['pro_measured'], y=result_unify['flux'], fit_reg=False)
plt.xlabel("Measured_protein_copy/cell")
plt.ylabel("Predicted_protein_copy/cell")







# part 2
# classify metabolic genes based on organelles
# then check the correlation in each organelles
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
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
organelle_v0 = ['mitochondrion', 'nucleus', 'cytosol',
 'endoplasmic reticulum', 'lipid droplet', 'fungal-type vacuole',
 'peroxisome', 'Golgi apparatus']
organelle_m0 = ['fungal-type vacuole membrane',
 'plasma membrane',
 'mitochondrial outer membrane',
 'endoplasmic reticulum membrane',
 'mitochondrial inner membrane',
 'Golgi membrane']
compartment_in0 = organelle_v0 + organelle_m0

# Note: it shows that four genes from 'Golgi membrane' were not related to the core metabolic functions from the model, so the predicted protein abudance is zero.


for organelle_target in compartment_in0:
    print(organelle_target)
    #test
    organelle_target = 'endoplasmic reticulum membrane'
    gene_target = m_gene_in_organelle[organelle_target]
    result_unify_c = result_unify[result_unify["geneID"].isin(gene_target)]
    plt.figure()
    sns.regplot(x=np.log10(result_unify_c['pro_measured'] + 1), y=np.log10(result_unify_c['flux'] + 1), fit_reg=False)
    plt.xlim(-0.5, 7)
    plt.ylim(-0.5, 7)
    plt.xlabel("log10(Measured_protein_copy/cell + 1)")
    plt.ylabel("log10(Predicted_protein_copy/cell +1)")
    plt.suptitle(organelle_target)
    # calculate the correlation coefficients - method2
    corr, ss = pearsonr(np.log10(result_unify_c['pro_measured'] + 1), np.log10(result_unify_c['flux'] + 1))
    print("Correlation coefficient:", corr)
    print("Correlation p_value:", ss)


# Note: It also shows that the predicted abundance of YJL167W is much higher than the measured one!
# 'endoplasmic reticulum' contains YJL167W
getRxnByGene(ecYeast, "YJL167W")
(sum(result_unify_c['flux'])-449824)/(sum(result_unify_c['pro_measured'])-84034)
(sum(result_unify_c['flux']))/(sum(result_unify_c['pro_measured']))

# Note: endoplasmic reticulum membrane, predicted protein abundance for YNR016C and YGR060W is much higher than measured