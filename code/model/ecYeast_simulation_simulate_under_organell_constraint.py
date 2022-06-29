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



# generate the general formula as the constraint
# all metabolic genes from ecGEMs
organelle_v = collectOrganelleTerm(type="volume")
organelle_m = collectOrganelleTerm(type="m")
gene_metabolic = gene_prot["geneID"].tolist()
compartment_in = organelle_v + organelle_m
m_gene_in_organelle = FingGenesForOrganelle(gene_set=gene_metabolic, compartment_list=compartment_in, compartment_type="organelle")


def SimulateOrganelleProAbundance(constraint_organelle, min_pro_abs, max_pro_abs, flux_expression, ecModel, saturation_cof = 0.44):
    # simulation in loop procedure
    # in vivo saturation of all, saturation_cof=0.44 for CENPK.113-7D strain
    lower = min_pro_abs * saturation_cof
    upper = max_pro_abs * saturation_cof
    step = (upper - lower) / 10
    organelle_abundance = []
    growth_rate = []
    predicted_usage = []
    ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)  # this is like the total protein pools in model, this constraint will affect growth prediction greatly.
    for upper_v in np.arange(lower, upper, step).tolist():
        model_tmp = ecModel.copy()
        same_flux = model_tmp.problem.Constraint(eval(flux_expression), lb=lower, ub=upper_v, name='same_flux')
        model_tmp.add_cons_vars(same_flux)
        # model_tmp = ecYeast # the model can't be used in the assignment
        # maximization
        objective = model_tmp.problem.Objective(
            model_tmp.reactions.r_4041.flux_expression,
            direction='max')  # biomass
        model_tmp.objective = objective
        solution2 = model_tmp.optimize()
        print("Max growth:" + str(solution2.objective_value))
        fluxes_select = solution2.fluxes[rxn_select]
        abundance_select0 = sum(list(fluxes_select))
        print("Total abundance:" + str(abundance_select0))
        organelle_abundance.append(upper_v)
        growth_rate.append(solution2.objective_value)
        predicted_usage.append(abundance_select0)

    # plot
    result_df = pd.DataFrame(
        {"growth": growth_rate, constraint_organelle: organelle_abundance, "predicted_usage": predicted_usage})
    plt.figure(figsize=(4, 4))
    sns.lineplot(x=constraint_organelle, y="growth", data=result_df, marker="o")
    plt.ylim(0, 0.5)
    plt.figure(figsize=(4, 4))
    sns.lineplot(x=constraint_organelle, y="predicted_usage", data=result_df, marker="o")

    # calculate the correlation coefficients
    from scipy.stats import pearsonr
    corr1, ss1 = pearsonr(organelle_abundance, growth_rate)
    corr2, ss2 = pearsonr(organelle_abundance, predicted_usage)
    return corr1, corr2, result_df


# loop for different organelle
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
 'peroxisomal membrane',
 'nuclear membrane']
compartment_in0 = organelle_v0 + organelle_m0
correlation1 = []
correlation2 = []
for org in compartment_in0:
    print(org)
    org = 'endoplasmic reticulum'
    compartment_info = organelle_pro_range[org].tolist()
    min_value = compartment_info[3]
    max_value = compartment_info[7]
    organelle_target = org
    gene_target = m_gene_in_organelle[organelle_target]
    rxn_select = gene_prot[gene_prot["geneID"].isin(gene_target)]["rxnID"].tolist()
    rxn_select = [x.replace("-A", "_A") for x in rxn_select]
    formula_list = ["model_tmp.reactions." + x + ".flux_expression" for x in rxn_select]
    formula_one = " + ".join(formula_list)
    c1, c2, df = SimulateOrganelleProAbundance(constraint_organelle=organelle_target, min_pro_abs=min_value, max_pro_abs=max_value, flux_expression=formula_one, ecModel=ecYeast)
    correlation1.append(c1)
    correlation2.append(c2)

result_df = pd.DataFrame({"compartment_in0": compartment_in0, "abundance_growth_cor": correlation1, "abundance_cor": correlation2})
result_df.to_excel("result/correlation_analysis.xlsx")


# check how saturation factor affect the model prediction
saturation_growth = {}
for saturation_cof0 in np.arange(0.2, 0.8, 0.1).tolist():
    org = 'endoplasmic reticulum'
    compartment_info = organelle_pro_range[org].tolist()
    min_value = compartment_info[3]
    max_value = compartment_info[7]
    organelle_target = org
    gene_target = m_gene_in_organelle[organelle_target]
    rxn_select = gene_prot[gene_prot["geneID"].isin(gene_target)]["rxnID"].tolist()
    rxn_select = [x.replace("-A", "_A") for x in rxn_select]
    formula_list = ["model_tmp.reactions." + x + ".flux_expression" for x in rxn_select]
    formula_one = " + ".join(formula_list)
    c1, c2, df = SimulateOrganelleProAbundance(constraint_organelle=organelle_target, min_pro_abs=min_value,
                                               max_pro_abs=max_value, flux_expression=formula_one, ecModel=ecYeast,
                                               saturation_cof=saturation_cof0)
    saturation_growth[saturation_cof0] = df['growth']


# analyze the result
df_saturation = pd.DataFrame(saturation_growth)
df_saturation.columns = ["ratio:" + str(round(x,2)) for x in np.arange(0.2, 0.8, 0.1)]
df_saturation['endoplasmic reticulum'] = df['endoplasmic reticulum']
# plot
# plot
plt.figure()
sns.lineplot(x='endoplasmic reticulum', y='value', hue='variable', style="variable",
             data=pd.melt(df_saturation, ['endoplasmic reticulum']))
plt.xlabel("endoplasmic reticulum's protein abundance (mmol/gDW)")
plt.ylabel("growth (/h)")







