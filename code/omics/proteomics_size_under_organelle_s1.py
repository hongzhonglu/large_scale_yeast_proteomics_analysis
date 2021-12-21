# this script is to process proteins with unit of molecular/cell directly recorded in literature.
# Hongzhong Lu
# 2021-11-20

import statistics

# import self function
from src.mainFunction import *
from src.protein_process import *

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

# input the pro compartment annotation data
pro_location = pd.read_excel("result/gene_compartment_mapping.xlsx")
pro_location = pro_location[['DBID', 'Systematic_name','GO_Name', 'Annot_Type', 'compartment']]

# input the pro abundance from two reference different sources
# the calculated result is quite similar to each other


# should further explore how to scale protein abundance to get absolute concentrations


# input data from SGD
pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')
pro_abundance = pro_abundance[["Systematic_name","Abundance_median", 'Abundance_SD']]
pro_abundance.columns = ["gene", "molecular/cell",'SD'] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]




# input data from cell system, 2018
pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
pro_abundance = pro_abundance[["Systematic Name","Mean molecules per cell","Median molecules per cell"]]
pro_abundance.columns = ["gene", "absolute_abundance","median_absolute_abundance"] # protein abundance per cell
pro_abundance.columns = ['gene','molecular/cell', "median_absolute_abundance"]
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]





# input data from nature methods 2014
pro_abundance = pd.read_excel("data/proteomics/data_Nature_method_2014.xlsx")
pro_abundance = pro_abundance[["ORF","Copy number"]]
pro_abundance.columns = ["gene", "absolute_abundance"] # protein abundance per cell
pro_abundance.columns = ['gene','molecular/cell']
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
pro_abundance = pro_abundance[pro_abundance["gene"].notna()]
pro_abundance = splitAbundance(pro_df=pro_abundance)




# test
# volume
location0 = 'mitochondrion'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
Vm,y= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)

location0 = 'nucleus'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
Vn,y= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)

location0 = 'cytosol'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
Vc,y= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)




# area
location0 = 'mitochondrial outer membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
x,Smom= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)

location0 = 'mitochondrial inner membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
x,Smim= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)

location0 = 'nuclear membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
x,Snm= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)

location0 = 'plasma membrane'
# for all related genes
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
x,Scellm= getStructureSize(pro_size0=pro_size, abundance0=pro_abundance1)