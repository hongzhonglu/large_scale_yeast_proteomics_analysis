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

    # using the updated version in 2024
    # Input the datasets from paxDB
    compartment = pd.read_csv("data/yeastmine_results_2024-04-15T10-39-56.tsv", sep='\t')
    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'Ontology_Description',
                           'GO_Namespace', 'GO_Name', 'GO_Identifier', 'Annot_Type', 'GO_Qualifier']
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

    # if filter == "Yes":
    #    return compartment_dict20
    # else:
    #    return compartment_dict_all0

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
    gene_plasma_membrane = pd.read_excel("data/sce_compartment_curation/plasma_membrane_annotations_v2.xlsx")
    gene_cell_wall = pd.read_excel("data/sce_compartment_curation/fungal_type_cell_wall_annotations_v3.xlsx")
    gene_fungal_type_vacuole_membrane = pd.read_excel(
        "data/sce_compartment_curation/fungal_type_vacuole_membrane_annotations_v2.xlsx")
    gene_nucleolus = pd.read_excel("data/sce_compartment_curation/nucleolus_annotations_v2.xlsx")
    gene_cytoplasm = pd.read_excel("data/sce_compartment_curation/cytoplasm_annotations_v2.xlsx")
    gene_cytosol = pd.read_excel("data/sce_compartment_curation/cytosol_annotations_v2.xlsx")
    gene_nucleus = pd.read_excel("data/sce_compartment_curation/nucleus_annotations_v2.xlsx")

    # mitochondrion specific
    gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx")
    # gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
    # gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/mitochondrion_xia.xlsx")
    gene_m_Outer_membrane = pd.read_excel(
        "data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v3.xlsx")
    gene_m_Inner_membrane = pd.read_excel(
        "data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v3.xlsx")
    gene_m_OI_space = pd.read_excel(
        "data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v3.xlsx")
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
            genes_select = [x for x in genes_select if
                            x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
        elif y == "endosome":
            genes_select = organelle1[y]
            genes_select = [x for x in genes_select if x not in [
                "YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
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

