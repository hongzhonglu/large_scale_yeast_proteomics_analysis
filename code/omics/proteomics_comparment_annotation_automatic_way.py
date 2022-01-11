# this script is to process compartment datasets
# 2021-11-16

# here the compartment annotation is mainly from SGD and MitoMiner
import pandas as pd

def getCompartmentGeneList(filter="Yes"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')

    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier',
                           'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
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
    annotation_type = list(set(annotation_type))
    # here if we remove "computational"
    compartment_with_evidence = compartment2[compartment2["Annot_Type"] != 'computational']
    compartment_with_computation = compartment2[compartment2["Annot_Type"] == 'computational']
    # in one procedure, if a protein has no compartment annotation from manual and high-throughput, then the computational is used!
    compartment_addition = compartment_with_computation[~compartment_with_computation["Systematic_name"].isin(compartment_with_evidence["Systematic_name"])]
    compartment_combine = pd.concat([compartment_with_evidence, compartment_addition])

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

    compartment_dict2 = {}
    for i, x in compartment_combine.iterrows():
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
    if filter == "Yes":
        return compartment_dict20
    else:
        return compartment_dict_all0



# test
compartment_dict_all0 = getCompartmentGeneList(filter="No")
compartment_dict20 = getCompartmentGeneList(filter="Yes")


# compare the difference
key0 = []
len0 = []

for key in compartment_dict_all0.keys():
    key0.append(key)
    len0.append(len(compartment_dict_all0[key]))
df1 = pd.DataFrame({"compartment":key0, "gene_number": len0})

key0 = []
len0 = []
for key in compartment_dict20.keys():
    key0.append(key)
    len0.append(len(compartment_dict20[key]))
df2 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# combine the dataframe
compartment_compare = pd.merge(left=df1, right=df2, left_on=['compartment'], right_on=['compartment'], how='left')
compartment_compare.columns = ["compartment", "all_annotation", "annotation_filter"]
compartment_compare.to_excel("data/compare_compartment_anotation_with_and_without_filter.xlsx")