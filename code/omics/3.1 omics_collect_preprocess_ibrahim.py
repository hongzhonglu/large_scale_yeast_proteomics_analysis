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
batch_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/batch_ExperimentalData.csv")
batch1 = batch_ExperimentalData.iloc[0,].tolist()
batch2 = list(batch_ExperimentalData.columns)
batch3 = [str(x) + "_miu=" + str(y)for x, y in zip(batch2, batch1)]
batch3[0] = "gene"
batch_ExperimentalData.columns = batch3
batch_ExperimentalData0 = batch_ExperimentalData.iloc[1:,]


translation_inhibit_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/TI_ExperimentalData.csv")
translation_inhibit1 = translation_inhibit_ExperimentalData.iloc[0,].tolist()
translation_inhibit2 = list(translation_inhibit_ExperimentalData.columns)
translation_inhibit3 = [str(x) + "_miu=" + str(y)for x, y in zip(translation_inhibit2, translation_inhibit1)]
translation_inhibit3[0] = "gene"
translation_inhibit_ExperimentalData.columns = translation_inhibit3
translation_inhibit_ExperimentalData0 = translation_inhibit_ExperimentalData.iloc[1:,]


chemostat_inhibit_ExperimentalData = pd.read_csv("data/ProteomicsData_ibrahim/chemostat_ExperimentalData.csv")
chemostat_inhibit1 = chemostat_inhibit_ExperimentalData.iloc[0,].tolist()
chemostat_inhibit2 = list(chemostat_inhibit_ExperimentalData.columns)
chemostat_inhibit3 = [str(x) + "_miu=" + str(y)for x, y in zip(chemostat_inhibit2, chemostat_inhibit1)]
chemostat_inhibit3[0] = "gene"
chemostat_inhibit_ExperimentalData.columns = chemostat_inhibit3
chemostat_inhibit_ExperimentalData0 = chemostat_inhibit_ExperimentalData.iloc[1:,]



# combine all the above dataset together
mass_fraction_ibrahim = pd.merge(left=batch_ExperimentalData0, right=translation_inhibit_ExperimentalData0, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_ibrahim = pd.merge(left=mass_fraction_ibrahim, right=chemostat_inhibit_ExperimentalData0, left_on=['gene'], right_on=['gene'], how="outer")



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
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_ibrahim, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_ibrahim.xlsx")





















