# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

import sys

# import self function
from src.mainFunction import *
from src.protein_process import *


# some general datasets
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

# absolute part
proteomics_NCB1 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_IO.xlsx", sheet_name="Table 10c. abs_prot_IO_SD108")

proteomics_NCB1 = proteomics_NCB1[["Entry","geneID","mean"]]
proteomics_NCB1.columns = ["Entry","geneID","IO_SD108_batch"]

# 1 relative data, change it as the absolute value
proteomics_NCB3 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_IO.xlsx", sheet_name="Table 11b. rel_prot_IO") # sce_FY4
proteomics_NCB4 = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_IO.xlsx", sheet_name="Table 11d. rel_prot_IO_resp_inh") # sce_cenpk
proteomics_NCB3 = proteomics_NCB3.iloc[:,0:17] # sce_FY4
proteomics_NCB4 = proteomics_NCB4.iloc[:,0:10] # sce_cenpk

proteomics_NCB3["ref_abs"] = multiMapping(proteomics_NCB1["IO_SD108_batch"],proteomics_NCB1["Entry"],proteomics_NCB3["Entry"])
proteomics_NCB3 = proteomics_NCB3.replace(to_replace='None', value=np.nan).dropna()
proteomics_NCB3['ref_abs'] = proteomics_NCB3['ref_abs'].apply(lambda x: float(x))
target_column_FY4 = list(proteomics_NCB3.columns)[2:17]
proteomics_NCB30 = proteomics_NCB3[["Entry"]+target_column_FY4]
for xx in target_column_FY4:
    print(xx)
    proteomics_NCB30[xx] = proteomics_NCB3.ref_abs * proteomics_NCB3[xx]
# unify the name
target_column_FY40 = ["IO_SD108_" + x for x in target_column_FY4]
miu_measured = ["miu=0.12", "miu=0.18", "miu=0.23", "miu=0.34","miu=0.45"]*3
target_column_FY40 = [x+"_"+y for x, y in zip(target_column_FY40, miu_measured)]


proteomics_NCB30.columns = ["Entry"] + target_column_FY40
# merge all the dataset
df_combine1 = pd.merge(left=proteomics_NCB1, right=proteomics_NCB30, left_on=['Entry'], right_on=['Entry'], how="left")
df_combine1 = df_combine1.rename(columns={'IO_SD108_batch': 'IO_SD108_batch_miu=0.52'})


# 2 relative data, change it as the absolute value
proteomics_NCB4["ref_abs"] = multiMapping(proteomics_NCB1["IO_SD108_batch"],proteomics_NCB1["Entry"],proteomics_NCB4["Entry"])
proteomics_NCB4 = proteomics_NCB4.replace(to_replace='None', value=np.nan).dropna()
proteomics_NCB4['ref_abs'] = proteomics_NCB4['ref_abs'].apply(lambda x: float(x))
target_column_CEN = list(proteomics_NCB4.columns)[1:10]
proteomics_NCB40 = proteomics_NCB4[["Entry"]+target_column_CEN]
for xx in target_column_CEN:
    print(xx)
    proteomics_NCB40[xx] = proteomics_NCB4.ref_abs * proteomics_NCB4[xx]
# unify the name
target_column_CEN0 = ["IO_SD108_" + x for x in target_column_CEN]
proteomics_NCB40.columns = ["Entry"] + target_column_CEN0
# merge the dataset
df_combine2 = pd.merge(left=proteomics_NCB1, right=proteomics_NCB40, left_on=['Entry'], right_on=['Entry'], how="left")
df_combine2 = df_combine2.rename(columns={'IO_SD108_batch': 'IO_SD108_batch_miu=0.52'})


# merge all the dataset
mass_fraction_NCB = pd.merge(left=df_combine1, right=df_combine2, left_on=['Entry'], right_on=['Entry'], how="outer")


mass_fraction_NCB = mass_fraction_NCB.drop(['geneID_x', 'geneID_y', 'IO_SD108_batch_miu=0.52_y'], axis=1)
mass_fraction_NCB = mass_fraction_NCB.rename(columns={'Entry': 'gene'})
mass_fraction_NCB = mass_fraction_NCB.rename(columns={'IO_SD108_batch_miu=0.52_x': 'IO_SD108_batch_miu=0.52'})

mass_fraction_final = mass_fraction_NCB.copy()
mass_fraction_final.to_excel("data/proteomics/mass_fraction_NCB_for_yeast_IO.xlsx")

def getCompartmentGeneList(filter="Yes"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    # compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t') # for sce
    compartment = pd.read_excel("data/nature_chemical_biology_datatset_2024/IO_gene_compartment.xlsx") # for yeast IO
    compartment = compartment.iloc[:,1:]

    # extract compartment
    compartment.columns = [ 'Systematic_name', 'GO_Name', 'Annot_Type']
    compartment["GO_Namespace"] = "cellular_component"
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
    # for compartment annotation removing some computation evidences
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

def ProMassRatio_Organelle(protein_abundance, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_abundance:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())

    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
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
            pro_abundance = protein_abundance[['gene', col0]]
            pro_abundance.columns = ['gene', 'g/gDW']
            '''
            # the following code is for sce specially!
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
            # get the sum
            pro_abundance.fillna(0, axis=1, inplace=True)
            pro_select = pro_abundance[pro_abundance['gene'].isin(genes_select)]
            sum_all = sum(pro_abundance['g/gDW'])
            sum_select = sum(pro_select['g/gDW'])
            ratio = sum_select/sum_all
            value1.append(ratio)
        result1[col0] = value1
    return result1



# test the above code
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_NCB_yeast_IO.xlsx")



