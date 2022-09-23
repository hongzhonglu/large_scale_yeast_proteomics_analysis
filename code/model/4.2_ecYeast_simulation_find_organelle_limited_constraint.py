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
def AddOrgConstraint(ecModel, flux_expression, min_pro_abs, max_pro_abs, constraint_name, saturation_cof = 0.44):
    model_tmp = ecModel.copy()
    lower = min_pro_abs * saturation_cof
    upper = max_pro_abs * saturation_cof
    model_tmp.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)
    same_flux = model_tmp .problem.Constraint(eval(flux_expression), lb=lower, ub=upper, name=constraint_name)
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
    #org = 'mitochondrial inner membrane' # just for the test
    compartment_info = organelle_pro_range[org].tolist()
    min_value = compartment_info[3] # minimum  value
    max_value = compartment_info[6] # 75%
    organelle_target = org
    gene_target = m_gene_in_organelle[organelle_target]
    rxn_select = gene_prot[gene_prot["geneID"].isin(gene_target)]["rxnID"].tolist()
    rxn_select = [x.replace("-A", "_A") for x in rxn_select]
    formula_list = ["model_tmp.reactions." + x + ".flux_expression" for x in rxn_select]
    formula_one = " + ".join(formula_list)
    constraint_name = org + '_constraint'
    constraint_name = constraint_name.replace(' ','_')
    ecYeast = AddOrgConstraint(ecModel=ecYeast, flux_expression=formula_one, min_pro_abs=min_value, max_pro_abs=max_value, constraint_name=constraint_name, saturation_cof = 0.44)
# check the growth
objective = ecYeast.problem.Objective(ecYeast.reactions.r_4041.flux_expression, direction='max') # biomass
ecYeast.objective = objective
solution2 = ecYeast.optimize()
print("Max growth:", solution2.objective_value)




# reset the constraints????
for org in compartment_in0:
    constraint_name = org + '_constraint'
    constraint_name = constraint_name.replace(' ','_')
    print(constraint_name)
    compartment_info = organelle_pro_range[org].tolist()
    min_value = compartment_info[3] # minimum  value
    max_value = compartment_info[6] # 75%
    ecYeast.constraints[constraint_name].ub = max_value
    ecYeast.constraints[constraint_name].lb = min_value
# check the growth
objective = ecYeast.problem.Objective(ecYeast.reactions.r_4041.flux_expression, direction='max') # biomass
ecYeast.objective = objective
solution2 = ecYeast.optimize()
print("Max growth:", solution2.objective_value)






# it found that if using the above constraint, the growth is very small. Some organelle protein total abundance is too strict.
# check the effect of constraints
ecYeast2 = ecYeast.copy() # copy model, each time only parameter is changed!
# reset the constraints???? Very strange that the constraint bounds changed when coping the models
for org in compartment_in0:
    constraint_name = org + '_constraint'
    constraint_name = constraint_name.replace(' ','_')
    print(constraint_name)
    compartment_info = organelle_pro_range[org].tolist()
    min_value = compartment_info[3] # minimum  value
    max_value = compartment_info[6] # 75%
    ecYeast2.constraints[constraint_name].ub = max_value
    ecYeast2.constraints[constraint_name].lb = min_value


org0 = 'endoplasmic reticulum membrane'
constraint_name0 = org0 + '_constraint'
constraint_name0 = constraint_name0.replace(' ','_')
print(constraint_name0)
compartment_info = organelle_pro_range[org0].tolist()
max_value = compartment_info[7]*4.5 # 9 for origninal ecYeast from DL
ecYeast2.constraints[constraint_name0].ub = max_value
objective = ecYeast2.problem.Objective(ecYeast2.reactions.r_4041.flux_expression, direction='max')  # biomass
ecYeast2.objective = objective
solution2 = ecYeast2.optimize()
flux_max = solution2.fluxes
print(solution2.objective_value)


# loops to find which constraint affect the model output
max_growth = []
for org in compartment_in0:
    with ecYeast2:
        constraint_name = org + '_constraint'
        constraint_name = constraint_name.replace(' ', '_')
        print(constraint_name)
        compartment_info = organelle_pro_range[org].tolist()
        max_value = compartment_info[6] * 2  # max value
        ecYeast2.constraints[constraint_name].ub = max_value
        # check which constraint affect the simulation??
        objective = ecYeast2.problem.Objective(ecYeast2.reactions.r_4041.flux_expression, direction='max')  # biomass
        ecYeast2.objective = objective
        solution2 = ecYeast2.optimize()
        max_growth.append(solution2.objective_value)

for x, y in zip(compartment_in0, max_growth):
    print(x, y)


