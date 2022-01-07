# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")

protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx")


# change the mmol/gDW as molecular/cell
coefficient1 = 7.8298e9
column0 = protein_abundance.columns
column1 = column0[1:]
protein_abundance1 = protein_abundance[column1]
protein_abundance2 = protein_abundance1.copy()
for x in column1:
    print(x)
    protein_abundance2[x] = protein_abundance1[x]*coefficient1

protein_abundance2['gene'] = protein_abundance['all_gene']


# combine the two dataframe
protein_copy_all = pd.merge(left=protein_abundance2, right=protein_copy, left_on=['gene'], right_on=['gene'], how='outer')
column2 = protein_copy_all.columns
column20 = [x for x in column2 if x is not 'gene']
column21 = ['gene'] + column20
protein_copy_all1 = protein_copy_all[column21]
protein_copy_all2 = protein_copy_all[column20]


ss = protein_copy_all2.describe()
protein_copy_all1.to_excel("data/proteomics/all_protein_copy.xlsx")
ss.to_excel("data/proteomics/protein_copy_statistical.xlsx")


# get the top 1000 proteins based on their molecular copies
new_df = protein_copy_all1[['gene']]
# Get the top 1000 proteins
for x in column20:
    print(x)
    df = protein_copy_all1[['gene',x]]
    ss2 = df.sort_values(by=[x], ascending=False)
    ss2_top1000 = ss2.iloc[0:1000, ]
    df[x][~df['gene'].isin(ss2_top1000['gene'])] = None
    new_df[x] = df[x]


# analysis
new_df1 = new_df[column20]
ss1 = new_df1.describe()
ss1.to_excel("data/proteomics/protein_copy_statistical_top1000.xlsx")



# get the size based on organelle and GO term

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

# input the pro compartment annotation data
pro_location = pd.read_excel("result/gene_compartment_mapping.xlsx")
pro_location = pro_location[['DBID', 'Systematic_name','GO_Name', 'Annot_Type', 'compartment']]


# area
organelle = ['mitochondrial outer membrane', 'mitochondrial inner membrane', 'nuclear membrane', 'plasma membrane']
result = pd.DataFrame({"compartment":organelle})
for col0 in column20:
    print(col0)
    value0=[]
    for y in organelle:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
        pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
        x, S = getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)
        value0.append(S)
    result[col0] = value0

result.to_excel("data/proteomics/membrance_size_across_conditions.xlsx")













# application 1 - check the glucose transporter abundance
glucose_transporter = pd.read_excel("data/glucose_transporter_abundance_check.xlsx")
glucose_transporter = glucose_transporter[glucose_transporter['Note'].isna()]
omics_glucose_transporter = protein_abundance[protein_abundance['all_gene'].isin(glucose_transporter['gene'])]
omics_glucose_transporter.to_excel("data/omics_glucose_transporter.xlsx")


glucose_transporter_copy = protein_copy[protein_copy['gene'].isin(glucose_transporter['gene'])]
glucose_transporter_copy['section_area'] = singleMapping(glucose_transporter['section_area'],glucose_transporter['gene'],glucose_transporter['gene'])
# get the section area
glucose_transporter_copy.to_excel("data/glucose_transporter_copy.xlsx")






