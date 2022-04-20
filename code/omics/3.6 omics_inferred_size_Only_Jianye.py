# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
column20 = list(protein_copy_all1.columns)
column20 = [x for x in column20 if x != "gene"]


# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]


# only analyze the jianye datasets under C limitation?
all_columns = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in all_columns if "_M" in x]
protein_copy_all1 = protein_copy_all1[Sample_ID_select+["gene"]]
compartment = getCompartmentGeneList(filter="Yes") # based on the automatic way
all_compartment = list(compartment.keys())
# compartment curation
# here for the following two membrane, only the compartment annotation with experimental evidence is used.
# also for the fungal-type vacuole membrane, some proteins belong to the metabolic enzymes with high abundance were removed.
# all_compartment = ['plasma membrane']
gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
# all_compartment = ['fungal-type vacuole membrane']
gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
result1 = pd.DataFrame({"compartment": all_compartment})
result2 = pd.DataFrame({"compartment": all_compartment})
for col0 in Sample_ID_select:
    print(col0)
    value1=[]
    value2=[]
    for y in all_compartment:
        print(y)
        location0 = y
        pro_abundance = protein_copy_all1[['gene',col0]]
        pro_abundance.columns = ['gene','molecular/cell']
        if y == "plasma membrane":
            genes_select = gene_plasma_membrane["gene"].tolist()# for the test
        elif y == "fungal-type vacuole membrane":
            genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
            genes_select = [x for x in genes_select if x not in ["YAL005C","YLL024C"]] # remove two genes for fungal type vacuole membrane
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
result1.to_excel("data/proteomics/volume_size_across_compartment_jianye_C_limitation.xlsx")
result2.to_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")