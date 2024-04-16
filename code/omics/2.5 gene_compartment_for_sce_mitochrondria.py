# note:
# this part can be rewritten as a function

import os
# import self function
from src.protein_process import *

# Part 1
# Initially check how many genes could find compartment
compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
#compartment = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
all_compartment = list(compartment.keys())
mitochondrion_related = ['mitochondrion',
'mitochondrial outer membrane',
'mitochondrial intermembrane space',
'mitochondrial matrix',
'mitochondrial inner membrane']


# get the subset
compartment_m = {x:y for x, y in compartment.items() if x in mitochondrion_related}
# based on the above annotation only a part of gene was assigned into each sub-organell of mitochondrion

# reanalyze the data set
mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx")
mitochondrion_gene = mitochondrion['gene'].tolist()

# combine
mitochondrial_inner_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_outer_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_matrix = pd.read_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v2.xlsx")

# combine
mitochondrial_intermembrane_space_annotations = pd.read_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v2.xlsx")




# input the NC paper
gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
gene_mitochondrion = gene_mitochondrion.iloc[0:987:]
gene_m_Outer_membrane = gene_mitochondrion[gene_mitochondrion['Outer membrane'] == "X"]
gene_m_Inner_membrane = gene_mitochondrion[gene_mitochondrion['Inner membrane'] == "X"]
gene_m_OI_space = gene_mitochondrion[gene_mitochondrion['Inter-membrane space'] == "X"]
gene_m_matrix = gene_mitochondrion[gene_mitochondrion['Matrix'] == "X"]



gene_m_Outer_membrane = list(set(gene_m_Outer_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_Inner_membrane = list(set(gene_m_Inner_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_OI_space = list(set(gene_m_OI_space['gene'].tolist()) & set(mitochondrion_gene))
gene_m_matrix = list(set(gene_m_matrix['gene'].tolist()) & set(mitochondrion_gene))



# further combine different data together:
gene_m_Outer_membrane_v3 = list(set(gene_m_Outer_membrane + mitochondrial_outer_membrane['gene'].tolist()))
gene_m_Inner_membrane_v3 = list(set(gene_m_Inner_membrane + mitochondrial_inner_membrane['gene'].tolist()))
gene_m_OI_space_v3 = list(set(gene_m_OI_space + mitochondrial_intermembrane_space_annotations['gene'].tolist()))
gene_m_matrix_v3 = list(set(gene_m_matrix + mitochondrial_matrix['gene'].tolist()))
gene_assigned = list(set(gene_m_Outer_membrane_v3 + gene_m_Inner_membrane_v3 + gene_m_OI_space_v3 + gene_m_matrix_v3)) # 742 proteins
gene_unassigned = list(set(mitochondrion_gene)-set(gene_assigned)) # 454 proteins



gene_m_Outer_membrane_v3 = pd.DataFrame({"gene": gene_m_Outer_membrane_v3})
gene_m_Inner_membrane_v3 = pd.DataFrame({"gene": gene_m_Inner_membrane_v3})
gene_m_OI_space_v3 = pd.DataFrame({"gene": gene_m_OI_space_v3})
gene_m_matrix_v3 = pd.DataFrame({"gene": gene_m_matrix_v3})
gene_m_unassigned_v3 = pd.DataFrame({"gene": gene_unassigned})


gene_m_Outer_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v3.xlsx")
gene_m_Inner_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v3.xlsx")
gene_m_OI_space_v3.to_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v3.xlsx")
gene_m_matrix_v3.to_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v3.xlsx")
gene_m_unassigned_v3.to_excel("data/sce_compartment_curation/mitochondrial_unassigned_manual_v3.xlsx")

