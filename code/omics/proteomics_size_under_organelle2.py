# this script is to try combine compartment annotation with protein reference abundance and size
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

# input the pro abundance from different sources
# should further explore how to scale protein abundance to get absolute concentrations
# input data from SGD
pro_abundance = pd.read_excel("data/yeast_proteomics_example_scale.xlsx")
pro_abundance = pro_abundance[["genes","copy_per_cell"]]
pro_abundance.columns = ["gene", "absolute_abundance"] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["absolute_abundance"].notna()]


# it shows after the data transformation, the calculate size of protein itself is larger than the organelle. So issues exist.
# first remove the proteins with abundance larger than 1000 0000 moleculars/cell.
abundance_cut_off = 8*1e7
pro_abundance = pro_abundance[pro_abundance['absolute_abundance'] <= abundance_cut_off]


# test
# volume
location0 = 'mitochondrion'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
Vm,y= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)

location0 = 'nucleus'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
Vn,y= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)

location0 = 'cytosol'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
Vc,y= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)





# area
location0 = 'mitochondrial outer membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
x,Smom= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)

location0 = 'mitochondrial inner membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
x,Smim= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)


location0 = 'nuclear membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
x,Snm= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)


location0 = 'plasma membrane'
genes_select = getGeneListFromLocation(gene_location_annotation=pro_location, location=location0)
x,Scellm= getStructureSize(genes_select0=genes_select,pro_size0=pro_size, pro_abundance0=pro_abundance)








# The following are two detailed examples
# with the mitochondrion volume as an example
location0 = 'mitochondrion'
genes_select = getGeneListFromLocation(pro_location, location0)
# change it as a dataframe
combine_df = pd.DataFrame({"gene": genes_select})
combine_df["Volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],combine_df["gene"])
combine_df["section_area"] = singleMapping(pro_size['section_area_new'],pro_size['locus'],combine_df["gene"])
combine_df["abundance"] = singleMapping(pro_abundance["absolute_abundance"],pro_abundance['gene'],combine_df["gene"])

# for the protein without abundance, use the median value from this group.
# calculate the abundance median value
abundance0 = combine_df["abundance"].tolist()
abundance1 = [x for x in abundance0 if np.isnan(x) == False]
abundance_median = statistics.median(abundance1)
abundance_update = []
for x in abundance0:
    if np.isnan(x) == False:
        x0 = x
    else:
        x0 = abundance_median
    abundance_update.append(x0)
combine_df["abundance_update"] = abundance_update

# calculate the size of all proteins in mitochondrion
# 1 纳米(nm)=0.001 微米(um)
total_volume = sum(combine_df["abundance_update"]*combine_df["Volume"])
# change nm^3 into um^3
total_volume_um = total_volume/1e9



# next we explore the area of sectional area calculated from proteins
# with the mitochondrion outer membrane as an example
# it initially shows that the organell membrane area is more constrainted than organell size
location0 = 'mitochondrial outer membrane'
genes_select = getGeneListFromLocation(pro_location, location0)
# change it as a dataframe
combine_df = pd.DataFrame({"gene": genes_select})
combine_df["Volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],combine_df["gene"])
combine_df["section_area"] = singleMapping(pro_size['section_area_new'],pro_size['locus'],combine_df["gene"])
combine_df["abundance"] = singleMapping(pro_abundance["absolute_abundance"],pro_abundance['gene'],combine_df["gene"])

# for the protein without abundance, use the median value from this group.
# calculate the abundance median value
abundance0 = combine_df["abundance"].tolist()
abundance1 = [x for x in abundance0 if np.isnan(x) == False]
abundance_median = statistics.median(abundance1)
abundance_update = []
for x in abundance0:
    if np.isnan(x) == False:
        x0 = x
    else:
        x0 = abundance_median
    abundance_update.append(x0)
combine_df["abundance_update"] = abundance_update

# calculate the size of all proteins in mitochondrion
# 1 纳米(nm)=0.001 微米(um)
total_area = sum(combine_df["abundance_update"]*combine_df["section_area"])
# change nm^2 into um^2
total_area_um = total_area/1e6