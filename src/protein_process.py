import math
import numpy as np
import statistics
import pandas as pd
from src.mainFunction import *
import os

def getSurfaceRatio(volume_ratio = 0.005/100, Vcell = 82, Scell=91.27):
    """
    The function is used to calculate the surface area ratio between organelle and cell.
    :param volume_ratio:
    :param Vcell:
    :param Scell:
    :return:
    """
    Vany = Vcell * volume_ratio
    Dany = 2 * (3 * Vany / (4 * math.pi)) ** (1 / 3)
    Sanym = 4 * math.pi * (Dany / 2) ** 2
    ratio = Sanym / Scell
    return ratio


def setReferenceProCopy():
    """
    This function is used to set the reference protein molecular/cell! If we still could not find the values from the
    reference value, then we can use the values at 5 percentile or 10 percentile?
    :return:
    """
    pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
    pro_abundance = pro_abundance[["Systematic Name", "Mean molecules per cell", "Median molecules per cell"]]
    pro_abundance.columns = ["gene", "absolute_abundance", "median_absolute_abundance"]  # protein abundance per cell
    pro_abundance.columns = ['gene', 'molecular/cell', "median_absolute_abundance"]
    pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
    reference_copy = pro_abundance
    # statistics_analysis = reference_copy.describe()
    v_five_percent = reference_copy['molecular/cell'].quantile(0.05)
    v_ten_percent = reference_copy['molecular/cell'].quantile(0.1)
    return reference_copy, v_five_percent, v_ten_percent


def getProAundance_old_version(genes_select0, pro_abundance0):
    """
    Note: this function need double check!!!

    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_abundance0: the unite is moleculars per cell
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    combine_df["molecular/cell"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],combine_df["gene"])

    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    combine_df.fillna(value=pd.np.nan, inplace=True) # change none into nan in the dataframe
    abundance0 = combine_df["molecular/cell"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    if len(abundance1) < 1:
        return "no_abundance"
    else:
        # use the first choice: for gene with no measured abundance, use the median value for the gene from the same compartment
        abundance_median = statistics.median(abundance1)  # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
        abundance_update = []
        for x in abundance0:
            if np.isnan(x) == False:
                x0 = x
            else:
                x0 = abundance_median
            abundance_update.append(x0)
        combine_df["molecular/cell_local"] = abundance_update


        # use the second choice: for gene with no measured abundance, use the value from the reference conditions??
        # load the reference molecular copies
        ref_abundance, v_5, v_10 = setReferenceProCopy()
        # set a dict
        gene_abundance = {}
        for i, x in ref_abundance.iterrows():
            print(i)
            gene_abundance[x['gene']] = x['molecular/cell']

        abundance_update2 = []
        for i, x in combine_df.iterrows():
            print(i)
            ss = x['molecular/cell']
            if np.isnan(ss) == False:
                x0 = ss
            elif x['gene'] in ref_abundance['gene'].tolist():
                x0 = gene_abundance[x['gene']]
            else:
                x0 = v_5
            abundance_update2.append(x0)

        combine_df["molecular/cell_global"] = abundance_update2

    return combine_df


def getProAundance(genes_select0, pro_abundance0):
    """
    Note: this function need double check!!!

    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_abundance0: the unite is moleculars per cell
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    #combine_df["molecular/cell"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],combine_df["gene"])
    combine_df = pd.merge(left=pro_abundance0, right=combine_df, left_on=['gene'], right_on=['gene'],how="right")


    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    combine_df.fillna(value=pd.np.nan, inplace=True) # change none into nan in the dataframe
    abundance0 = combine_df["molecular/cell"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    if len(abundance1) < 1:
        return "no_abundance"
    else:

        """
        # use the first choice: for gene with no measured abundance, use the median value for the gene from the same compartment
        abundance_median = statistics.median(abundance1)  # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
        abundance_update = []
        for x in abundance0:
            if np.isnan(x) == False:
                x0 = x
            else:
                x0 = abundance_median
            abundance_update.append(x0)
        combine_df["molecular/cell_local"] = abundance_update


        # use the second choice: for gene with no measured abundance, use the value from the reference conditions??
        # load the reference molecular copies
        ref_abundance, v_5, v_10 = setReferenceProCopy()
        # set a dict
        gene_abundance = {}
        for i, x in ref_abundance.iterrows():
            print(i)
            gene_abundance[x['gene']] = x['molecular/cell']

        abundance_update2 = []
        for i, x in combine_df.iterrows():
            print(i)
            ss = x['molecular/cell']
            if np.isnan(ss) == False:
                x0 = ss
            elif x['gene'] in ref_abundance['gene'].tolist():
                x0 = gene_abundance[x['gene']]
            else:
                x0 = v_5
            abundance_update2.append(x0)

        combine_df["molecular/cell_global"] = abundance_update2

    return combine_df"""
        return combine_df


def getStructureSize(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0
    #combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    #combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])

    # the following three lines were used to update the upper one
    combine_df = pd.merge(combine_df, pro_size0, left_on="gene", right_on="locus", how="left")
    combine_df.columns = ['gene', 'molecular/cell', 'DBID', 'locus', 'Volume', 'section_area']
    combine_df = combine_df[['gene', 'molecular/cell', 'Volume', 'section_area']]

    #combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])
    # it shows that some genes have no locus
    combine_df = combine_df[~combine_df["section_area"].isna()]

    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell_global"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell_global"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_volume_um, total_area_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell_global"] * combine_df["Volume"]
        combine_df["total_area"] = combine_df["molecular/cell_global"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_volume_um, total_area_um, combine_df

def getStructureSize_MeasuredAbundances(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0

    #combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    #combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])

    # the following three lines were used to update the upper two
    combine_df = pd.merge(combine_df, pro_size0, left_on="gene", right_on="locus", how="left")
    combine_df.columns = ['gene', 'molecular/cell', 'DBID', 'locus', 'Volume', 'section_area']
    combine_df = combine_df[['gene', 'molecular/cell', 'Volume', 'section_area']]

    #combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])
    # it shows that some genes have no locus
    combine_df = combine_df[~combine_df["section_area"].isna()]
    combine_df = combine_df[~combine_df["molecular/cell"].isna()] # newly added for this new function


    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_volume_um, total_area_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell"] * combine_df["Volume"]
        combine_df["total_area"] = combine_df["molecular/cell"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_volume_um, total_area_um, combine_df

def get_total_protein_volume(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """
    # should make sure no structure size data is nan
    combine_df = abundance0
    #combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])

    # the following three lines were used to update the upper one
    combine_df = pd.merge(combine_df,pro_size0, left_on="gene", right_on="locus", how="left")
    combine_df.columns = ['gene', 'molecular/cell', 'DBID', 'locus', 'Volume','section_area']
    combine_df = combine_df[['gene', 'molecular/cell', 'Volume','section_area']]

    combine_df = combine_df[~combine_df["molecular/cell"].isna()] # newly added for this new function
    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    if need_check=="No":
        return total_volume_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell"] * combine_df["Volume"]
        return total_volume_um, combine_df




# get the compartments of all genes

# now we have the updated version of the compartment
def getCompartmentGeneList(type="all"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    #compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')
    # extract compartment
    #compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier', 'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']


    # using the updated version in 2024
    # Input the datasets from paxDB
    compartment = pd.read_csv("data/yeastmine_results_2024-04-15T10-39-56.tsv", sep='\t')
    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'Ontology_Description', 'GO_Namespace', 'GO_Name', 'GO_Identifier', 'Annot_Type', 'GO_Qualifier']
    compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]

    # filter out compartment with "complex" or "subunit"
    compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("subunit")]
    # firstly remove some general cellular component
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("cellular_component")]

    # analyze the annotation type
    annotation_type = compartment2["Annot_Type"].tolist()
    print("Annotation type:")
    print(list(set(annotation_type)))


    # here if we remove "computational"
    # it seems that this method is wrong, as a protein could be located in multiple compartment. Thus it removes too much information
    compartment_with_evidence = compartment2[compartment2["Annot_Type"] != 'computational']
    compartment_with_computation = compartment2[compartment2["Annot_Type"] == 'computational']

    len(set(compartment_with_evidence["Systematic_name"].tolist()))
    len(set(compartment_with_computation["Systematic_name"].tolist()))

    # in one procedure, if a protein has no compartment annotation from manual and high-throughput, then the computational is used!
    #compartment_addition = compartment_with_computation[~compartment_with_computation["Systematic_name"].isin(compartment_with_evidence["Systematic_name"])]
    #compartment_combine = pd.concat([compartment_with_evidence, compartment_addition])


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


    # for compartment annotation from manual evidence
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


    #if filter == "Yes":
    #    return compartment_dict20
    #else:
    #    return compartment_dict_all0
    
    if type == "all":
        return compartment_dict_all0
    elif type =="manual":
        return compartment_dict20
    else:
        return compartment_dict30


# calibrate the compartments of some genes based on manual experiment
# however this step can be omitted.
def getCompartment_manual_curation():
    # design a function to process the original compartment annotation from SGD
    data_dir = "/Users/xluhon/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/sce_compartment_curation/original_annotation/"
    protein_complex = pd.read_excel("data/complex_info.xlsx")
    protein_complex["subunit"] = protein_complex["subunit"].str.replace("-MONOMER", "")
    all_file0 = os.listdir(data_dir)
    for xx in all_file0:
        print(xx)
        try:
            ss0 = open(data_dir + xx).readlines()
            index0 = [i for i, x in enumerate(ss0) if "Gene/Complex" in x]
            ss0_update = ss0[index0[0] + 1:]
            # get the gene list
            gene_list = [x.split("\t")[1] for x in ss0_update]
            # classify it as two types
            list1 = [x for x in gene_list if "CPX-" in x]
            list2 = list(set([x for x in gene_list if "CPX-" not in x]))
            # get the ID from the complex
            complex_select = protein_complex[protein_complex['complex'].isin(list1)]
            if complex_select.shape[0] > 1:
                gene_belong_complex = complex_select['subunit'].tolist()
                list_new = list(set(list2 + gene_belong_complex))
            else:
                list_new = list2

            new_df = pd.DataFrame({"gene": list_new})
            new_df.to_excel(
                "/Users/xluhon/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/sce_compartment_curation/" + xx.replace(
                    ".txt", "_v2.xlsx"))
        except:
            pass
#getCompartment_manual_curation()

# when running the following function, please firstly run function - getCompartment_manual_curation
gene_cell_wall = pd.read_excel("data/sce_compartment_curation/fungal_type_cell_wall_annotations_v2.xlsx")
gene_anotation = pd.read_excel("data/sce_gene_annotation_SGD.xlsx")
gene_anotation_for_cell_wall = gene_anotation[gene_anotation['ID'].isin(gene_cell_wall['gene'])]
# here 9 enzyme genes belong to cell wall were removed.
#gene_anotation_for_cell_wall.to_excel("data/sce_compartment_curation/fungal_type_cell_wall_annotations_v3.xlsx")
# getCompartment_manual_curation() # run the compartment curation preprocess. If update the compartment information, need to run this function

def gene_location_curation_sce(organelle0):
    # use some manually checked gene compartment definion
    # if the manual curated gene number for one compartment is larger, nealy equal to computational, then use the manual curation
    # otherwise using the computation prediction???
    # input the annotation from sgd
    # organelle0 = getCompartmentGeneList(type="all") # this is just for the test
    gene_plasma_membrane = pd.read_excel("data/sce_compartment_curation/plasma_membrane_annotations_v2.xlsx")
    gene_cell_wall = pd.read_excel("data/sce_compartment_curation/fungal_type_cell_wall_annotations_v3.xlsx")
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/sce_compartment_curation/fungal_type_vacuole_membrane_annotations_v2.xlsx")
    gene_nucleolus = pd.read_excel("data/sce_compartment_curation/nucleolus_annotations_v2.xlsx")
    gene_cytoplasm = pd.read_excel("data/sce_compartment_curation/cytoplasm_annotations_v2.xlsx")
    gene_cytosol = pd.read_excel("data/sce_compartment_curation/cytosol_annotations_v2.xlsx")
    gene_nucleus = pd.read_excel("data/sce_compartment_curation/nucleus_annotations_v2.xlsx")

    # mitochondrion specific
    gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx")
    #gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
    #gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrion_xia.xlsx")
    gene_m_Outer_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v3.xlsx")
    gene_m_Inner_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v3.xlsx")
    gene_m_OI_space = pd.read_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v3.xlsx")
    gene_m_matrix = pd.read_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v3.xlsx")
    gene_m_unassigned = pd.read_excel("data/sce_compartment_curation/mitochondrial_unassigned_manual_v3.xlsx")

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
        else:
            genes_select = organelle1[y]
        organelle0_update[y] = list(filter(lambda x: str(x) != 'nan', genes_select))
        # remove cytoplasm
        # organelle0_update = {x:y for x, y in organelle0_update.items() if x is not "cytoplasm"}
    return organelle0_update


def get_total_membrane_area(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0

    # combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])

    # the following three lines were used to update the upper one
    combine_df = pd.merge(combine_df, pro_size0, left_on="gene", right_on="locus", how="left")
    combine_df.columns = ['gene', 'molecular/cell', 'DBID', 'locus', 'Volume', 'section_area']
    combine_df = combine_df[['gene', 'molecular/cell', 'Volume', 'section_area']]

    # it shows that some genes have no locus
    combine_df = combine_df[~combine_df["section_area"].isna()]
    combine_df = combine_df[~combine_df["molecular/cell"].isna()] # newly added for this new function
    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_area_um
    else:
        combine_df["total_area"] = combine_df["molecular/cell"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_area_um, combine_df

# calculate the mass ratio
def ProMassRatio_Organelle(protein_abundance, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_abundance:
    :param compartment_type:
    :return:
    """
    # test
    # protein_abundance = mass_fraction_final

    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment) # based on the SGD manual curation
        all_compartment = list(compartment.keys())

    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    # use some manually checked gene compartment definion
    # gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    # gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    # creat a dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            # test
            # y = "plasma membrane"
            # col0 = "Glucose_phase_rep1(g/gDW)"
            pro_abundance = protein_abundance[['gene', col0]]
            pro_abundance.columns = ['gene', 'g/gDW']

            '''if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
            genes_select = compartment[y]

            # get the sum
            pro_abundance.fillna(0, axis=1, inplace=True)
            pro_select = pro_abundance[pro_abundance['gene'].isin(genes_select)]
            sum_all = sum(pro_abundance['g/gDW'])
            sum_select = sum(pro_select['g/gDW'])
            ratio = sum_select/sum_all
            value1.append(ratio)
        result1[col0] = value1
    return result1

# calculate the volume ratio
def Pro_3D_Volume_Ratio_Cal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment) # based on the SGD manual curation
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    #gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    #gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")

    # creat two dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    #result2 = pd.DataFrame({"compartment": all_compartment})

    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        pro_abundance = protein_copy[['gene', col0]]
        pro_abundance.columns = ['gene', 'molecular/cell']
        total_volume = get_total_protein_volume(pro_size0=pro_size, abundance0=pro_abundance, need_check="No")

        for y in all_compartment:
            print(y)
            # test
            # y = "cytosol"
            
            '''
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
            genes_select = compartment[y]


            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
            if pro_abundance1 is "no_abundance":
                value1.append(None)
            else:
                x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
                value1.append(x/total_volume)
                #value2.append(S)
        result1[col0] = value1
    return result1

# check in the following two functions, whether the cell wall proteins are covered!
def getMembraneProList():
    # for the membrane annotation, uniprot is good in transmembrane annotation
    # maybe get the intersection between uniprot and SGD in membrane annotation
    # also, put the unassigned membrane protein as a unique group?

    # be careful about this part of analysis. There are so many membrane proteins without manual curation!!
    # as the first step: define the membrane or transporter protein
    #protein_transporter = open("/Users/xluhon/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/tcdb.txt").readlines()
    #protein_transporter = [x for x in protein_transporter if ">" in x]
    #protein_transporter = [x for x in protein_transporter if "S288c" in x]
    #protein_ID = []
    #for xx in protein_transporter:
    #    ss0 = xx.split("|")[2]
    #    protein_ID.append(ss0)
    # get the transporter gene id in sce
    #uniprotGeneID_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
    #transporter_tf = uniprotGeneID_mapping[uniprotGeneID_mapping["Entry"].isin(protein_ID)]
    #transporter_pro_list = transporter_tf["GeneName"].tolist()

    # get the membrane annotation from SGD
    membrane_pro = pd.read_excel("data/sce_compartment_curation/membrane_annotations_computational_v2.xlsx")
    membrane_pro_list = list(set(membrane_pro['gene'].tolist()))

    # Input the datasets from paxDB
    #compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')
    # extract compartment
    #compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier',
    #                       'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
    #compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]
    # using the updated version in 2024
    # Input the datasets from paxDB
    compartment = pd.read_csv("data/yeastmine_results_2024-04-15T10-39-56.tsv", sep='\t')
    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'Ontology_Description', 'GO_Namespace', 'GO_Name', 'GO_Identifier', 'Annot_Type', 'GO_Qualifier']
    compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]
    compartment1_membrane_filter = compartment1[compartment1["GO_Name"].str.contains("membrane")][compartment1["GO_Name"] !="mitochondrial intermembrane space"]
    membrane_pro_list_database = list(set(compartment1_membrane_filter["Systematic_name"].tolist()))

    # check the relation between the annotation from the above procedures
    membrane_pro_final_merge = list(set(membrane_pro_list) & set(membrane_pro_list_database))
    cell_wall = pd.read_excel("data/sce_compartment_curation/fungal_type_cell_wall_annotations_v2.xlsx")
    cell_wall_gene = cell_wall['gene'].tolist()
    # plus transporter proteins and cell wall proteins
    # but the total membrane proteins are too many!!???
    # membrane_pro_final_merge11 = list(set(transporter_pro_list) - set(membrane_pro_final_merge)) + membrane_pro_final_merge + list(set(cell_wall_gene) - set(membrane_pro_final_merge))
    membrane_pro_final_merge11 = membrane_pro_final_merge + list(set(cell_wall_gene) - set(membrane_pro_final_merge))

    return list(set(membrane_pro_final_merge11))

# # calculate the membrane ratio
def Pro_Membrance_Ratio_Cal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment) # based on the SGD manual curation
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    #gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    #gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")

    # creat two dataframe to save the result
    all_compartment = [x for x in all_compartment if "membrane" in x]
    all_compartment = [x for x in all_compartment if x != "mitochondrial intermembrane space"] + ['fungal-type cell wall']
    result2 = pd.DataFrame({"compartment": all_compartment})
    membrane_pro_final_merge = getMembraneProList()
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        value2 = []
        # for the total membrane
        pro_abundance = protein_copy[['gene', col0]]
        pro_abundance.columns = ['gene', 'molecular/cell']
        all_membrane_abundance1 = getProAundance(genes_select0=membrane_pro_final_merge, pro_abundance0=pro_abundance)
        S_total = get_total_membrane_area(pro_size0=pro_size, abundance0=all_membrane_abundance1)

        for y in all_compartment:
            print(y)
            '''
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
            genes_select = compartment[y]


            genes_select = list(set(genes_select) & set(membrane_pro_final_merge))
            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)

            if pro_abundance1 is "no_abundance":
                value1.append(None)
                value2.append(None)
            else:
                x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
                value2.append(S/S_total)
        result2[col0] = value2
    return result2


# absolute protein structure volume and area for each organelle
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
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment) # based on the SGD manual curation
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    #gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    #gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")

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
            '''
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
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

# GO-term level calculation
def Pro3DCalForGOterm(protein_copy):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]

    go_term = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")
    all_go_term = list(go_term.keys())
    result1 = pd.DataFrame({"go_term": all_go_term})
    result2 = pd.DataFrame({"go_term": all_go_term})
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        value2 = []
        for y in all_go_term:
            print(y)
            # location0 = y
            pro_abundance = protein_copy[['gene', col0]]
            pro_abundance.columns = ['gene', 'molecular/cell']
            genes_select = go_term[y]
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

# absolute protein abundance for each organelle
def ProAbsoluteCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment) # based on the SGD manual curation
        all_compartment = list(compartment.keys())
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    #gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    #gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
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
            '''
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
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

# absolute protein abundance for each go_term
def ProAbsoluteCal2(protein_copy, compartment_type="go_term"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "go_term":
        # go term
        go_term = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")
        all_go_term = list(go_term.keys())
        # sample ID information
        Sample_ID_select = list(protein_copy.columns)
        Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
        out = pd.DataFrame({"go_term": all_go_term})
        # result2 = pd.DataFrame({"go_term": all_go_term})
        for col0 in Sample_ID_select:
            print(col0)
            value1 = []
            # value2=[]
            for y in all_go_term:
                print(y)
                # location0 = y
                pro_abundance = protein_copy[['gene', col0]]
                pro_abundance.columns = ['gene', 'molecular/cell']
                genes_select = go_term[y]
                pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
                if pro_abundance1 is "no_abundance":
                    value1.append(None)
                    # value2.append(None)
                else:
                    pro_abundance1 = pro_abundance1.dropna()
                    x = sum(pro_abundance1['molecular/cell'])
                    value1.append(x)
                    # value2.append(S)
            out[col0] = value1
        return out


def splitAbundance(pro_df):
    """
    The function is used to quality check of protein abundance in molecular/cell or other unit before entering next step.
    :param pro_df: A dataframe should columns-gene,molecular/cell.
    :return:
    """
    # sometimes it shows mutiple proteins together have one abundance value, here we need a function to do the quality check!
    len1 = pro_df.shape[0]
    colnames=pro_df.columns
    pro_df1 = pro_df[pro_df['gene'].str.contains(';')]
    if (len(pro_df1) > 0):
        print('Mutiple protein have one abudance value! Need quality check.')
    pro_df2 = pro_df[~pro_df['gene'].str.contains(';')]
    gene0 = []
    abundance0 = []
    for i, x in pro_df1.iterrows():
        print(i, x)
        s = x["gene"].split(";")
        len0 = len(s)
        gene0 = gene0 + s
        v = [x[colnames[1]] / len0] * len0
        abundance0 = abundance0 + v
    gene0 = [x.strip(" ") for x in gene0]
    pro_df1 = pd.DataFrame({"gene": gene0, colnames[1]: abundance0})

    pro_df = pd.concat([pro_df1, pro_df2], axis=0)
    len2 = pro_df.shape[0]

    if (len2 > len1):
        print('Complete the quality check!')

    return pro_df


def getGoTermGeneList(input1, input2):
    """
    This function is generate the GO term and its gene list
    :parameter1: a directory contain a excel file with gene GO term annotation
    :parameter2: a directory contain a excel file with gene id mapping
    :return: A dict. GO_term as key and gene list as value
    ------------
    Usage:
    GO_term_gene = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")

    Hongzhong Lu
    2022.01.07
    """
    # Input the datasets from paxDB
    GO_term = pd.read_excel(input1)
    GO_term['GO-slim mapper process term'] = GO_term['GO-slim mapper process term'].str.strip()

    gene_info = pd.read_csv(input2, sep="\t")
    gene_short_name = gene_info['gene_name'].tolist()
    gene_short_name = [str(x) for x in gene_short_name]
    gene_info['gene_name'] = gene_short_name
    gene_short_name2 = []
    gene_locus = gene_info['locus'].tolist()
    for x, y in zip(gene_short_name, gene_locus):
        print(x, y)
        if x == 'nan':
            gene_short_name2.append(y)
        else:
            gene_short_name2.append(x)
    # build a dict
    # it shows that some genes have no locus
    # be careful in this step
    gene_name_dict = {}
    for w, v in zip(gene_short_name2, gene_locus):
        gene_name_dict[w] = v

    # build GO_term dict
    GO_dict = {}
    for i, x in GO_term.iterrows():
        print(i)
        ss = list(x)
        name = ss[0]
        ss = ss[1:]
        mylist = [str(x) for x in ss]
        newlist = [v for v in mylist if v != 'nan']
        # get the gene OFR name base on short gene ID
        key_all = gene_name_dict.keys()
        newlist1 = []
        for x in newlist:
            if x in key_all:
                x0 = gene_name_dict[x]
            else:
                x0 = x
            newlist1.append(x0)
        GO_dict[name] = newlist1
    return GO_dict



def getGeneListFromLocation(gene_location_annotation, location):
    """
    The function is to extract gene list based on its compartment information
    :param gene_location_annotation: A dataframe contains the annotation of each protein
    :param location: A string represent the compartment name
    :return:
    """
    #location = 'mitochondrial envelope'
    gene_subset = gene_location_annotation[gene_location_annotation["GO_Name"]==location]
    gene_list = list(set(gene_subset["Systematic_name"].tolist()))
    return gene_list




def AllProteomicsAnalysis(pro_df):
    """
    The function is used to do the general statistical analysis of proteomics datasets across conditions.
    :param pro_df: A dataframe to store proteomics, with column "gene"
    :return:

    _____

    usage: AllProteomicsAnalysis(pro_df=protein_copy_all1)
    """
    protein_copy_all2 = pro_df.drop(columns=['gene'])
    ss = protein_copy_all2.describe()
    # get the top 1000 proteins based on their molecular copies
    new_df = pro_df[['gene']]
    # Get the top 1000 proteins
    column20 = protein_copy_all2.columns
    for x in column20:
        print(x)
        df = pro_df[['gene', x]]
        ss2 = df.sort_values(by=[x], ascending=False)
        ss2_top1000 = ss2.iloc[0:1000, ]
        df[x][~df['gene'].isin(ss2_top1000['gene'])] = None
        new_df[x] = df[x]
    # analysis
    new_df1 = new_df[column20]
    ss1 = new_df1.describe()
    ss1.to_excel("data/proteomics/protein_copy_statistical_top1000.xlsx")
    ss.to_excel("data/proteomics/protein_copy_statistical.xlsx")
    return ss


def calculateCoefficient(cell_volume0):
    cell_volume = cell_volume0  # 32.6 # fL/cell
    dry_content = 0.35  # https://onlinelibrary.wiley.com/doi/pdf/10.1002/j.2050-0416.1952.tb02660.x#:~:text=Yeast%20cakes%20produced%20by%20normal,the%20conditions%20of%20growth%2C%20to
    cell_density = 1.1126e-12  # g/fL yeast cell density under exponential growth, [g/fL] = 1e12  g/mL
    coefficent10 = 1000 / 6.022e+23 / cell_volume / dry_content / cell_density  # from molecular/cell into mmol/gDW
    coefficent20 = 1 / coefficent10  # from mmol/gDW into molecular/cell
    return coefficent20


def calculateCurationCoefficent():
    # curation of rosemary datasets based on the fitted cell volume under different growth rates
    # fitting formula to calculate the coefficients
    # when miu > 0.2, cell_volume = 77.32 miu + 15.771
    # when miu < 0.2, average volume is 28 um^3

    # Rosemary sample ID
    Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
    coefficient2 = 7.8298e9  # this is the original coefficient used in cell system paper!
    growth_rate = [0.05, 0.1, 0.13, 0.18, 0.3, 0.35]
    #cell_volume_fit = [28, 28, 28, 28, 38.967, 42.833] # here assume there exist a minimum cell size
    cell_volume_fit = [x*47.458+22.742 for x in growth_rate] # here we assume there exist a linear increase of cell size when growth rate increased
    coefficent_list = [calculateCoefficient(cell_volume0=x) for x in cell_volume_fit]
    curation_coefficent = [x / coefficient2 for x in coefficent_list]
    new_coefficent = []
    for x in curation_coefficent:
        print(x)
        s = [x] * 3
        new_coefficent = new_coefficent + s
    cell_volume_all = []
    for x in cell_volume_fit:
        print(x)
        s = [x] * 3
        cell_volume_all = cell_volume_all + s
    growth_all = []
    for x in growth_rate:
        print(x)
        s = [x] * 3
        growth_all = growth_all + s

    curation_info_rosemary = pd.DataFrame({"ID":Sample_ID_select,"growth_rate":growth_all, "cell_size":cell_volume_all, "curation_coefficent": new_coefficent})

    return curation_info_rosemary


def collectOrganelleTerm(type):
    """
    Some compartment need manual check.
    The function is just to get the important organelle list for volume or membrane size calculation.
    :param type:
    :return:

    usage:
    organelle_v = collectOrganelleTerm(type="volume")
    organelle_m = collectOrganelleTerm(type="m")

    """
    volume_list = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum', 'endosome', 'lipid droplet',
                   'fungal-type vacuole', 'peroxisome', 'ribosome', 'Golgi apparatus', 'cytosolic ribosome',
                   'mitochondrial ribosome', 'nucleolus']
    membrane_list = ['fungal-type vacuole membrane', 'plasma membrane', 'mitochondrial outer membrane',
                     'prospore membrane', 'endoplasmic reticulum membrane', 'mitochondrial inner membrane',
                     'Golgi membrane', 'cellular bud membrane',
                     'late endosome membrane', 'peroxisomal membrane', 'nuclear membrane', 'endosome membrane',
                     'nuclear inner membrane']
    if type == "volume":
        return volume_list
    else:
        return membrane_list


def FingGenesForOrganelle(gene_set, compartment_list, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        # all_compartment = list(compartment.keys())

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    all_compartment = compartment_list
    result_df = dict()
    for y in all_compartment:
            print(y)
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            # here we need calculate the intersection
            result_df[y] = list(set(genes_select) & set(gene_set))
    return result_df

'''
def Pro3DCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
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
    return result1, result2'''

'''
def ProAbsoluteCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
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
    return result1'''

