# note:
# this part can be rewritten as a function
import pandas as pd
import os
# import self function
from src.protein_process import *
# reanalyze the data set
def getCompartmentGeneList(type="all"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    # compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')
    # extract compartment
    # compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier', 'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
    # input gene annotation in Uniprot
    sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

    # using the updated version in 2026
    # Input the datasets from paxDB
    compartment = pd.read_csv("data/alliancemine_results_2026-01-03T10-33-37.tsv", sep='\t')
    compartment = compartment.iloc[:,0:11]
    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name',  'Standard_name', 'Feature_type', 'Gene_qualifier','GO_Identifier', 'GO_Name', 'GO_Namespace', 'Annot_Type', 'GO_Qualifier', 'Organism']

    compartment = compartment[compartment['Systematic_name'].isin(sce_gene['GeneName'])]
    compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]

    # filter out compartment with "complex" or "subunit"
    compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
    compartment2 = compartment2[~compartment1["DBID"].str.contains("CPX-")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("subunit")]
    # firstly remove some general cellular component
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("cellular_component")]

    # more filter
    compartment2 = compartment2[compartment2['GO_Qualifier']=='located_in']  # only select "located"


    # analyze the annotation type
    annotation_type = compartment2["Annot_Type"].tolist()
    print("Annotation type:")
    print(list(set(annotation_type)))
    ratios = pd.Series(annotation_type).value_counts(normalize=True).round(4)
    print(ratios)

    # add annotation
    go_evidence_codes = {
    'IDA': 'Inferred from Direct Assay',
    'IEA': 'Inferred from Electronic Annotation',
    'HDA': 'Inferred from High Throughput Direct Assay',
    'IBA': 'Inferred from Biological aspect of Ancestor',
    'TAS': 'Traceable Author Statement',
    'ISA': 'Inferred from Sequence Alignment',
    'ND': 'No biological Data available',  # no any evidence
    'IPI': 'Inferred from Physical Interaction',
    'NAS': 'Non-traceable Author Statement',
    'ISS': 'Inferred from Sequence or structural Similarity',
    'IMP': 'Inferred from Mutant Phenotype',
    'IC': 'Inferred by Curator',
    'EXP': 'Inferred from Experiment',
    'IGI': 'Inferred from Genetic Interaction',
    'ISM': 'Inferred from Sequence Model',
    'ISO': 'Inferred from Sequence Orthology'}

    exp = [
    'IDA',  # Inferred from Direct Assay (直接实验测定)
    'HDA',  # Inferred from High Throughput Direct Assay (高通量直接实验)
    'IPI',  # Inferred from Physical Interaction (物理相互作用实验)
    'IMP',  # Inferred from Mutant Phenotype (突变体表型实验)
    'IGI',  # Inferred from Genetic Interaction (遗传相互作用实验)
    'IC',   # Inferred by Curator (人工curator推断)
    'EXP',  # Inferred from Experiment (实验证据伞状代码)
    'TAS',  # Traceable Author Statement (可追溯作者声明，已弃用)
    'NAS',  # Non-traceable Author Statement (不可追溯作者声明，已弃用)
     ]


    comput = [
    'IEA',  # Inferred from Electronic Annotation (纯电子/自动注释)
    'IBA',  # Inferred from Biological aspect of Ancestor (系统发育祖先推断)
    'ISA',  # Inferred from Sequence Alignment (序列比对)
    'ISS',  # Inferred from Sequence or structural Similarity (序列/结构相似性，广义)
    'ISM',  # Inferred from Sequence Model (序列模型)
    'ISO'   # Inferred from Sequence Orthology (正交同源性)
    ]

    # here if we remove "computational"
    # it seems that this method is wrong, as a protein could be located in multiple compartment. Thus it removes too much information
    compartment_with_evidence = compartment2[compartment2["Annot_Type"].isin(exp)]
    compartment_with_computation = compartment2[compartment2["Annot_Type"].isin(comput)]

    len(set(compartment_with_evidence["Systematic_name"].tolist()))
    len(set(compartment_with_computation["Systematic_name"].tolist()))

    # in one procedure, if a protein has no compartment annotation from manual and high-throughput, then the computational is used!
    # compartment_addition = compartment_with_computation[~compartment_with_computation["Systematic_name"].isin(compartment_with_evidence["Systematic_name"])]
    # compartment_combine = pd.concat([compartment_with_evidence, compartment_addition])

    # just for the test
    # compartment_combine_test = compartment_combine[compartment_combine['GO_Name'] == "nuclear membrane"]

    # build the dict
    compartment_dict_all = {}
    for i, x in compartment2.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict_all.keys():
            compartment_dict_all[x['GO_Name']] = list(set(compartment_dict_all[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict_all[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict_all0 = {}
    for key in compartment_dict_all.keys():
        print(key)
        value = compartment_dict_all[key]
        if len(value) >= 6:
            compartment_dict_all0[key] = value
        else:
            pass

    # for compartment annotation from manual evidence
    compartment_dict2 = {}
    for i, x in compartment_with_evidence.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict2.keys():
            compartment_dict2[x['GO_Name']] = list(set(compartment_dict2[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict2[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict20 = {}
    for key in compartment_dict2.keys():
        print(key)
        value = compartment_dict2[key]
        if len(value) >= 6:
            compartment_dict20[key] = value
        else:
            pass

    # for compartment annotation from computational evidence
    compartment_dict3 = {}
    for i, x in compartment_with_computation.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict3.keys():
            compartment_dict3[x['GO_Name']] = list(set(compartment_dict3[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict3[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict30 = {}
    for key in compartment_dict3.keys():
        print(key)
        value = compartment_dict3[key]
        if len(value) >= 6:
            compartment_dict30[key] = value
        else:
            pass
    if type == "all":
        return compartment_dict_all0
    elif type == "manual":
        return compartment_dict20
    else:
        return compartment_dict30

def gene_location_curation_sce(organelle0):
    # use some manually checked gene compartment definion
    # if the manual curated gene number for one compartment is larger, nealy equal to computational, then use the manual curation
    # otherwise using the computation prediction???
    # input the annotation from sgd
    # organelle0 = getCompartmentGeneList(type="all") # this is just for the test
    gene_plasma_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/plasma_annotations.xlsx")
    gene_cell_wall = pd.read_excel("data/sce_compartment_curation/2026_curated/fungal_type_cell_wall_annotations.xlsx")
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/fungal_type_vacuole_membrane_annotations.xlsx")
    gene_nucleolus = pd.read_excel("data/sce_compartment_curation/2026_curated/nucleolus_annotations.xlsx")
    gene_cytoplasm = pd.read_excel("data/sce_compartment_curation/2026_curated/cytoplasm_annotations.xlsx")
    gene_cytosol = pd.read_excel("data/sce_compartment_curation/2026_curated/cytosol_annotations.xlsx")
    gene_nucleus = pd.read_excel("data/sce_compartment_curation/2026_curated/nucleus_annotations.xlsx")
    gene_peroxisome = pd.read_excel("data/sce_compartment_curation/2026_curated/peroxisome_annotations.xlsx")

    gene_endoplasmic_reticulum = pd.read_excel("data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_annotations.xlsx")
    gene_endoplasmic_reticulum_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_membrane_annotations.xlsx")
    gene_Golgi_apparatus = pd.read_excel("data/sce_compartment_curation/2026_curated/Golgi_apparatus_annotations.xlsx")
    gene_Golgi_membrane_annotations = pd.read_excel("data/sce_compartment_curation/2026_curated/Golgi_membrane_annotations.xlsx")





    # mitochondrion specific
    gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrion_annotations.xlsx")
    gene_m_Outer_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_outer_membrane_annotations.xlsx")
    gene_m_Inner_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_inner_membrane_annotations.xlsx")
    gene_m_OI_space = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_intermembrane_space_annotations.xlsx")
    gene_m_matrix = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_matrix_annotations.xlsx")
    gene_m_unassigned = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_unassigned.xlsx")

    organelle1 = organelle0.copy()
    organelle1['mitochondrion_unassigned'] = gene_m_unassigned['gene'].tolist()

    organelle0_update = {}
    for y in organelle1.keys():
        print(y)
        if y == "plasma membrane":
            genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
        elif y == "fungal-type cell wall":
            genes_select = gene_cell_wall["gene"].tolist()
        elif y == "fungal-type vacuole membrane":
            genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
            genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
        elif y == "endosome":
            genes_select = organelle1[y]
            genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
        elif y == "nucleolus":
            genes_select = gene_nucleolus["gene"].tolist()  # for the test
        elif y == "cytoplasm":
            genes_select = gene_cytoplasm["gene"].tolist()  # for the test
        elif y == "cytosol":
            genes_select = gene_cytosol["gene"].tolist()  # for the test
        elif y == "mitochondrion":
            genes_select = gene_mitochondrion["gene"].tolist()  # for the test
        elif y == "mitochondrial outer membrane":
            genes_select = gene_m_Outer_membrane["gene"].tolist()  # for the test
        elif y == "mitochondrial inner membrane":
            genes_select = gene_m_Inner_membrane["gene"].tolist()  # for the test
        elif y == "mitochondrial intermembrane space":
            genes_select = gene_m_OI_space["gene"].tolist()  # for the test
        elif y == "mitochondrial matrix":
            genes_select = gene_m_matrix["gene"].tolist()  # for the test
        elif y == "nucleus":
            genes_select = gene_nucleus["gene"].tolist()  # for the test
        elif y == "peroxisome":
            genes_select = gene_peroxisome["gene"].tolist()  # for the test
        elif y == "endoplasmic reticulum":
            genes_select = gene_endoplasmic_reticulum["gene"].tolist()  # for the test
        elif y == "endoplasmic reticulum membrane":
            genes_select = gene_endoplasmic_reticulum_membrane["gene"].tolist()  # for the test
        elif y == "Golgi apparatus":
            genes_select = gene_Golgi_apparatus["gene"].tolist()  # for the test
        elif y == "Golgi membrane":
            genes_select = gene_Golgi_membrane_annotations["gene"].tolist()  # for the test
        else:
            genes_select = organelle1[y]
        organelle0_update[y] = list(filter(lambda x: str(x) != 'nan', genes_select))
       # remove duplicates
        for key in organelle0_update:
            organelle0_update[key] = list(dict.fromkeys(organelle0_update[key]))

    return organelle0_update


# compartment info
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
compartment_corrected = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
# remove duplicates
compartment = {key: list(set(value)) for key, value in compartment.items()}
compartment_corrected = {key: list(set(value)) for key, value in compartment_corrected.items()}

# check cytosol
c1 = compartment['cytosol']
c2 = compartment_corrected['cytosol']
print(set(c1)-set(c2))
# note: cytosol中不少蛋白证据为TAS，但缺乏进一步文献支撑。

# compare the difference
key0 = []
len0 = []
for key in compartment.keys():
    key0.append(key)
    len0.append(len(compartment[key]))
df1 = pd.DataFrame({"compartment":key0, "gene_number": len0})


key0 = []
len0 = []
for key in compartment_corrected.keys():
    key0.append(key)
    len0.append(len(compartment_corrected[key]))
df2 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# combine the dataframe
compartment_compare = pd.merge(left=df1, right=df2, left_on=['compartment'], right_on=['compartment'], how='outer')
compartment_compare.columns = ["compartment", "annotation_combine", "annotation_curation"]
compartment_compare.to_excel("data/compare_compartment_annotation_with_and_without_manual_curation.xlsx")

# save the corrected compartment annotation
mapping =[]
for key, value in compartment_corrected.items():
    print(key, value)
    new0 = [key+"@"+ x for x in value]
    mapping = mapping + new0
df0 = pd.DataFrame({"pair": mapping})
df3 = df0['pair'].str.split('@', n=1, expand=True)
df3.columns = ['compartment','gene']
df3.to_excel("data/compartment_sce_curation.xlsx")





