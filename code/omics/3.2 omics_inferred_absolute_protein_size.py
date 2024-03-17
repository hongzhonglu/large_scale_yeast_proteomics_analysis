# Note:
# all these analysis is based on the protein copy?

# import self function
from src.model_process import *
from src.protein_process import *



# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")

# only analyze the rosemary datasets under NH4 limitation?
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all_rosemary = protein_copy_all1[Sample_ID_select+["gene"]]
# curate the protein copies based on newly calculated size
df_curated = calculateCurationCoefficent()
curation_coefficent = df_curated["curation_coefficent"].to_list()
for i, sid in enumerate(Sample_ID_select):
    print(i, sid, curation_coefficent[i])
    protein_copy_all_rosemary[sid] = protein_copy_all_rosemary[sid]*curation_coefficent[i]

protein_copy_all_rosemary.to_excel("data/proteomics/all_protein_copy_rosemary.xlsx")

s1, s2 = Pro3DCal(protein_copy_all_rosemary)
s1.to_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")


s1, s2 = Pro3DCalForGOterm(protein_copy_all_rosemary)
s1.to_excel("data/proteomics/volume_size_across_go_term_Rosemary.xlsx")
s2.to_excel("data/proteomics/membrance_size_across_go_term_Rosemary.xlsx")

# absolute protein abundance for each organelle or go-term
result1 = ProAbsoluteCal(protein_copy_all_rosemary,compartment_type="organelle")
result1.to_excel("data/proteomics/total_protein_abundance_across_compartment_Rosemary_NH4_limitation_v2.xlsx")

result2 = ProAbsoluteCal2(protein_copy=protein_copy_all_rosemary, compartment_type="go_term")
result2.to_excel("data/proteomics/total_protein_abundance_across_compartment_Rosemary_NH4_limitation_v2.xlsx")










# jianye datasets
all_columns = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in all_columns if "_M" in x]
protein_copy_jianye = protein_copy_all1[Sample_ID_select+["gene"]]
s1, s2 = Pro3DCal(protein_copy_jianye)
s1.to_excel("data/proteomics/volume_size_across_compartment_jianye_C_limitation.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")


# all datasets
s1, s2 = Pro3DCal(protein_copy_all1)
s1.to_excel("data/proteomics/volume_size_across_compartment.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")

s1, s2 = Pro3DCalForGOterm(protein_copy_all1)
s1.to_excel("data/proteomics/volume_size_across_go_term.xlsx")
s2.to_excel("data/proteomics/membrance_size_across_go_term.xlsx")


# Tao2 datasets
protein_copy_all1 = pd.read_excel("data/proteomics/protein_copy_tao_nc.xlsx")
s1, s2 = Pro3DCal(protein_copy_all1)
s1.to_excel("data/proteomics/volume_size_across_compartment_tao2.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment_tao2.xlsx")




## next we can focus on genes only from the ecGEMs to set the additional constraints
# also select other samples to calculate the general trend in enzyme volume from each organelle
# input the protein abundance data
# here it should be careful that Jianye datasets are not calibrated based on the cell size information!!
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
# based on all 1150 genes
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gene_list = []
for gene in ecYeast.genes:
    print(gene.id)
    gene_list.append(gene.id)
protein_copy_all_select = protein_copy_all1[protein_copy_all1["gene"].isin(gene_list)]

s1, s2 = Pro3DCal(protein_copy_all_select)
# for metabolic enzyme
s1.to_excel("data/proteomics/ecGEM_volume_size_across_compartment.xlsx")
s2.to_excel("data/proteomics/ecGEM_membrane_size_across_compartment.xlsx")





# calculate protein abundance for each organelle in unit of mmol/gDCW
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx") # the unit the g/gDW
# test
omics_combine_auto.columns = omics_combine_auto.columns.str.replace('all_gene', 'gene')

# based on all 1150 genes
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gene_list = []
for gene in ecYeast.genes:
    print(gene.id)
    gene_list.append(gene.id)
omics_combine_auto_select = omics_combine_auto[omics_combine_auto["gene"].isin(gene_list)]

s111 = ProAbsoluteCal(omics_combine_auto_select)
# for metabolic enzyme
s111.to_excel("data/proteomics/ecGEM_absolute_pro_across_compartment.xlsx")


