# Note:
# all these analysis is based on the protein copy?

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
                genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
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
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())
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
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
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



# all datasets
s1, s2 = Pro3DCal(protein_copy_all1)
s1.to_excel("data/proteomics/volume_size_across_compartment.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")



s1 = pd.read_excel("data/proteomics/volume_size_across_compartment.xlsx")
s2 = pd.read_excel("data/proteomics/membrane_size_across_compartment.xlsx")



# calculate the total volume of proteins
# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

protein_copy_all1["pro_volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],protein_copy_all1['gene'])
sample_ID = list(protein_copy_all1.columns)
sample_ID = sample_ID[1:77]
volume_list = []
for x in sample_ID:
    print(x)
    ss1 = protein_copy_all1[[x,"pro_volume"]]
    ss1["value"] = ss1[x]*ss1["pro_volume"]
    ss1 = ss1[~ss1["value"].isna()]
    sum0 = sum(ss1['value'])
    # change nm^3 into um^3
    total_volume_um = sum0 / 1e9
    volume_list.append(total_volume_um)
# creat a new dataframe
total_pro_volume = pd.DataFrame({"sampleID":sample_ID,"total_pro_volume":volume_list})

# calculate the fraction of organell protein volume per total volume

s1_matrix = s1.iloc[:,1:77]

s1_matrix_new = s1_matrix
for i in range(0,76,1):
    print(i)
    s1_matrix_new.iloc[:, i] = s1_matrix.iloc[:, i]/volume_list[i]

s1_matrix_new['compartment'] = s1['compartment']

s1_matrix_new.to_excel("data/proteomics/compartment_volume_fraction.xlsx")
