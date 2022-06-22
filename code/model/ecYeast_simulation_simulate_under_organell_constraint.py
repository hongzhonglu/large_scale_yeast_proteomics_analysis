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
same_flux = ecYeast.problem.Constraint(ecYeast.reactions.prot_YLR342W.flux_expression + ecYeast.reactions.prot_YPR159W.flux_expression + ecYeast.reactions.prot_YJL005W.flux_expression + ecYeast.reactions.prot_YCR024C_A.flux_expression + ecYeast.reactions.prot_YGL008C.flux_expression + ecYeast.reactions.prot_YPL036W.flux_expression + ecYeast.reactions.prot_YGR060W.flux_expression + ecYeast.reactions.prot_YNL192W.flux_expression + ecYeast.reactions.prot_YOR171C.flux_expression + ecYeast.reactions.prot_YOR317W.flux_expression + ecYeast.reactions.prot_YMR246W.flux_expression + ecYeast.reactions.prot_YJL100W.flux_expression + ecYeast.reactions.prot_YLR305C.flux_expression + ecYeast.reactions.prot_YDR208W.flux_expression + ecYeast.reactions.prot_YMR008C.flux_expression + ecYeast.reactions.prot_YLR020C.flux_expression + ecYeast.reactions.prot_YOR348C.flux_expression + ecYeast.reactions.prot_YCR010C.flux_expression + ecYeast.reactions.prot_YER056C.flux_expression + ecYeast.reactions.prot_YER060W.flux_expression + ecYeast.reactions.prot_YER060W_A.flux_expression + ecYeast.reactions.prot_YGL186C.flux_expression + ecYeast.reactions.prot_YJR152W.flux_expression + ecYeast.reactions.prot_YDR384C.flux_expression + ecYeast.reactions.prot_YGR121C.flux_expression + ecYeast.reactions.prot_YNL142W.flux_expression + ecYeast.reactions.prot_YGR065C.flux_expression + ecYeast.reactions.prot_YDR536W.flux_expression + ecYeast.reactions.prot_YLR081W.flux_expression + ecYeast.reactions.prot_YKL217W.flux_expression + ecYeast.reactions.prot_YDR342C.flux_expression + ecYeast.reactions.prot_YDR345C.flux_expression + ecYeast.reactions.prot_YHR092C.flux_expression + ecYeast.reactions.prot_YHR094C.flux_expression + ecYeast.reactions.prot_YMR011W.flux_expression + ecYeast.reactions.prot_YOR011W.flux_expression + ecYeast.reactions.prot_YCR098C.flux_expression + ecYeast.reactions.prot_YGL084C.flux_expression + ecYeast.reactions.prot_YLL043W.flux_expression + ecYeast.reactions.prot_YCL025C.flux_expression + ecYeast.reactions.prot_YKR039W.flux_expression + ecYeast.reactions.prot_YOL020W.flux_expression + ecYeast.reactions.prot_YDR497C.flux_expression + ecYeast.reactions.prot_YOL103W.flux_expression + ecYeast.reactions.prot_YMR319C.flux_expression + ecYeast.reactions.prot_YMR058W.flux_expression + ecYeast.reactions.prot_YBR068C.flux_expression + ecYeast.reactions.prot_YEL063C.flux_expression + ecYeast.reactions.prot_YBR069C.flux_expression + ecYeast.reactions.prot_YCR075C.flux_expression + ecYeast.reactions.prot_YGR191W.flux_expression + ecYeast.reactions.prot_YGR055W.flux_expression + ecYeast.reactions.prot_YCR028C.flux_expression + ecYeast.reactions.prot_YBR296C.flux_expression + ecYeast.reactions.prot_YCR037C.flux_expression + ecYeast.reactions.prot_YJL198W.flux_expression + ecYeast.reactions.prot_YJL129C.flux_expression + ecYeast.reactions.prot_YLL028W.flux_expression + ecYeast.reactions.prot_YOR273C.flux_expression + ecYeast.reactions.prot_YPL274W.flux_expression + ecYeast.reactions.prot_YLL061W.flux_expression + ecYeast.reactions.prot_YLR138W.flux_expression + ecYeast.reactions.prot_YHL016C.flux_expression + ecYeast.reactions.prot_YGR138C.flux_expression + ecYeast.reactions.prot_YPR156C.flux_expression + ecYeast.reactions.prot_YBR294W.flux_expression + ecYeast.reactions.prot_YLR092W.flux_expression + ecYeast.reactions.prot_YPL092W.flux_expression + ecYeast.reactions.prot_YLR237W.flux_expression + ecYeast.reactions.prot_YOR071C.flux_expression + ecYeast.reactions.prot_YBR021W.flux_expression + ecYeast.reactions.prot_YBL042C.flux_expression + ecYeast.reactions.prot_YLL052C.flux_expression + ecYeast.reactions.prot_YPR192W.flux_expression + ecYeast.reactions.prot_YNL065W.flux_expression + ecYeast.reactions.prot_YOR306C.flux_expression + ecYeast.reactions.prot_YMR162C.flux_expression + ecYeast.reactions.prot_YBR295W.flux_expression + ecYeast.reactions.prot_YML125C.flux_expression + ecYeast.reactions.prot_YDR038C.flux_expression + ecYeast.reactions.prot_YDR039C.flux_expression + ecYeast.reactions.prot_YDR040C.flux_expression + ecYeast.reactions.prot_YJR040W.flux_expression + ecYeast.reactions.prot_YKL220C.flux_expression + ecYeast.reactions.prot_YLR214W.flux_expression + ecYeast.reactions.prot_YNR060W.flux_expression + ecYeast.reactions.prot_YOL152W.flux_expression + ecYeast.reactions.prot_YOR381W.flux_expression + ecYeast.reactions.prot_YKR093W.flux_expression + ecYeast.reactions.prot_YDR093W.flux_expression + ecYeast.reactions.prot_YOL122C.flux_expression + ecYeast.reactions.prot_YNL275W.flux_expression + ecYeast.reactions.prot_YLR130C.flux_expression + ecYeast.reactions.prot_YOL130W.flux_expression,
    lb=0.00002,
    ub=0.00006)
ecYeast.add_cons_vars(same_flux)


same_flux = ecYeast.problem.Constraint(
    eval(formula_one),
    lb=0.00002,
    ub=0.00006)
ecYeast.add_cons_vars(same_flux)





