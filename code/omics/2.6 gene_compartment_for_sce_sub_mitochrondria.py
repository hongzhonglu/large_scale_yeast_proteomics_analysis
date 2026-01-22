# note:
# this part can be rewritten as a function

import os
import pandas as pd

# import self function
from src.protein_process import *

# input gene annotation in Uniprot
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")



# reanalyze the data set
gene_mitochondrion_sgd = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrion_annotations.xlsx") # this is based sgd annotation
mitochondrion_gene = list(set(gene_mitochondrion_sgd['gene'].tolist()))


# combine
def polish_annotaiton(df):
    df = df[df['Qualifier']=='located in'] # only select "located"
    df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
    # remove the protein with IEA evidence. IEA (Electronic Annotation)	最低	仅基于序列相似性预测，未经验证。
    df = df[~df['Evidence'].isin(['IEA'])]
    df_renamed = df.rename(columns={'Systematic Name/Complex Accession': 'gene'})
    return df_renamed

mitochondrial_inner_membrane = pd.read_csv("data/sce_compartment_curation/2026/mitochondrial_inner_membrane_annotations.txt",sep='\t', skiprows=8, header=0)
mitochondrial_inner_membrane = polish_annotaiton(mitochondrial_inner_membrane)

mitochondrial_outer_membrane = pd.read_csv("data/sce_compartment_curation/2026/mitochondrial_outer_membrane_annotations.txt",sep='\t', skiprows=8, header=0)
mitochondrial_outer_membrane = polish_annotaiton(mitochondrial_outer_membrane)

mitochondrial_matrix = pd.read_csv("data/sce_compartment_curation/2026/mitochondrial_matrix_annotations.txt",sep='\t', skiprows=8, header=0)
mitochondrial_matrix = polish_annotaiton(mitochondrial_matrix)

mitochondrial_intermembrane_space_annotations = pd.read_csv("data/sce_compartment_curation/2026/mitochondrial_intermembrane_space_annotations.txt",sep='\t', skiprows=8, header=0)
mitochondrial_intermembrane_space_annotations = polish_annotaiton(mitochondrial_intermembrane_space_annotations)




# input the NC paper
gene_mitochondrion_nc = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
gene_mitochondrion_nc = gene_mitochondrion_nc.iloc[0:987:]
gene_m_Outer_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Outer membrane'] == "X"]
gene_m_Inner_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Inner membrane'] == "X"]
gene_m_OI_space = gene_mitochondrion_nc[gene_mitochondrion_nc['Inter-membrane space'] == "X"]
gene_m_matrix = gene_mitochondrion_nc[gene_mitochondrion_nc['Matrix'] == "X"]


# 不同方法的交集
gene_m_Outer_membrane = list(set(gene_m_Outer_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_Inner_membrane = list(set(gene_m_Inner_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_OI_space = list(set(gene_m_OI_space['gene'].tolist()) & set(mitochondrion_gene))
gene_m_matrix = list(set(gene_m_matrix['gene'].tolist()) & set(mitochondrion_gene))



# 不同方法的并集
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


gene_m_Outer_membrane_v3.to_excel("data/sce_compartment_curation/2026_curated/mitochondrial_outer_membrane_annotations.xlsx")
gene_m_Inner_membrane_v3.to_excel("data/sce_compartment_curation/2026_curated/mitochondrial_inner_membrane_annotations.xlsx")
gene_m_OI_space_v3.to_excel("data/sce_compartment_curation/2026_curated/mitochondrial_intermembrane_space_annotations.xlsx")
gene_m_matrix_v3.to_excel("data/sce_compartment_curation/2026_curated/mitochondrial_matrix_annotations.xlsx")
gene_m_unassigned_v3.to_excel("data/sce_compartment_curation/2026_curated/mitochondrial_unassigned.xlsx")