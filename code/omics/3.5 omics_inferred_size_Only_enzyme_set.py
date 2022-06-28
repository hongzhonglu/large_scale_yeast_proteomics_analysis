# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle
# In other condition, the protein copy/cell under different conditions are calibrated based on the cell size datasets.


# import self function
from src.model_process import *
from src.protein_process import *

# TO-DO: the function will be further optimized so that it could be used other analysis
def Pro3DCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")

    # creat two dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    result2 = pd.DataFrame({"compartment": all_compartment})

    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        value2 = []
        for y in all_compartment:
            print(y)
            pro_abundance = protein_copy[['gene', col0]]
            pro_abundance.columns = ['gene', 'molecular/cell']
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
            if pro_abundance1 is "no_abundance":
                value1.append(None)
                value2.append(None)
            else:
                x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
                value1.append(x)
                value2.append(S)
        result1[col0] = value1
        result2[col0] = value2
    return result1, result2


def ProAbsoluteCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    # creat a dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            pro_abundance = protein_copy[['gene', col0]]
            pro_abundance.columns = ['gene', 'molecular/cell']
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
            if pro_abundance1 is "no_abundance":
                value1.append(None)
            else:
                pro_abundance1 = pro_abundance1.dropna()
                x = sum(pro_abundance1['molecular/cell'])
                value1.append(x)
        result1[col0] = value1
    return result1


# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
# only analyze the rosemary datasets under NH4 limitation?
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all1 = protein_copy_all1[Sample_ID_select+["gene"]]
protein_copy_all_rosemary = protein_copy_all1.copy()
# curate the protein copies based on newly calculated size
df_curated = calculateCurationCoefficent()
curation_coefficent = df_curated["curation_coefficent"].to_list()
for i, sid in enumerate(Sample_ID_select):
    print(i, sid, curation_coefficent[i])
    protein_copy_all_rosemary[sid] = protein_copy_all_rosemary[sid]*curation_coefficent[i]

protein_copy_all_rosemary.to_excel("data/proteomics/all_protein_copy_rosemary.xlsx")





"""
# read the model
GEM_yeast = read_sbml_model('/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xml')
gene_GEM = getALLGEMgene()
protein_copy_all_rosemary = protein_copy_all_rosemary[protein_copy_all_rosemary["gene"].isin(gene_GEM)]
# for metabolic enzyme
s1, s2 = Pro3DCal(protein_copy_all_rosemary)
s1.to_excel("data/proteomics/GEM_volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")
s2.to_excel("data/proteomics/GEM_membrane_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")
"""



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






## input the absolute protein abundance
# protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")
protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx")

# based on all 1150 genes
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gene_list = []
for gene in ecYeast.genes:
    print(gene.id)
    gene_list.append(gene.id)
protein_abundance_all_select = protein_abundance[protein_abundance["all_gene"].isin(gene_list)]
protein_abundance_all_select = protein_abundance_all_select.rename(columns={'all_gene': 'gene'})
p1 = ProAbsoluteCal(protein_abundance_all_select)
p1.to_excel("data/proteomics/ecGEM_absolute_pro_across_compartment.xlsx", index=False)
