#  this script is to process proteins with unit of molecular/cell directly calculated from literature.
# Hongzhong Lu
# 2021-11-20

import statistics

# import self function
from src.mainFunction import *
from src.protein_process import *
from src.model_process import *

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

# input the pro compartment annotation data
pro_location = pd.read_excel("result/gene_compartment_mapping.xlsx")
pro_location = pro_location[['DBID', 'Systematic_name','GO_Name', 'Annot_Type', 'compartment']]



# input the pro abundance from different sources

# input dataset2
pro_abundance = pd.read_excel("data/proteomics/data_PNAS_2021_scale.xlsx")
pro_abundance = pro_abundance[["gene","molecular/cell"]]
pro_abundance.columns = ['gene','molecular/cell']
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]



# it shows after the data transformation, the calculate size of protein itself is larger than the organelle. So issues exist.
# first remove the proteins with abundance larger than 1000 0000 moleculars/cell.
abundance_cut_off = 8*1e7
pro_abundance = pro_abundance[pro_abundance['molecular/cell'] <= abundance_cut_off]
pro_abundance = splitAbundance(pro_df=pro_abundance)


gene_GEM = getALLGEMgene()


# test
location0 = 'plasma membrane'
# for all related genes
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
x1, Scellm1, analysis_df1 = getStructureSize(pro_size0=pro_size, abundance0=abundance1, need_check="Yes")


# for the metabolic genes
genes_select_model = list(set(genes_select) & set(gene_GEM))
abundance2 = getProAundance(genes_select0=genes_select_model, pro_abundance0=pro_abundance)
x2, Scellm2, analysis_df2 = getStructureSize(pro_size0=pro_size, abundance0=abundance2, need_check="Yes")


genes_select_all_transporter = getTransporterCellMembraneGEM()
abundance3 = getProAundance(genes_select0=genes_select_all_transporter, pro_abundance0=pro_abundance)
x3, Scellm3, analysis_df3 = getStructureSize(pro_size0=pro_size, abundance0=abundance3, need_check="Yes")


genes_select_glucose = getProteinForRxnGEM(rxnID=['r_1166'])
abundance4 = getProAundance(genes_select0=genes_select_glucose, pro_abundance0=pro_abundance)
x4, Scellm4, analysis_df4 = getStructureSize(pro_size0=pro_size, abundance0=abundance4, need_check="Yes")





