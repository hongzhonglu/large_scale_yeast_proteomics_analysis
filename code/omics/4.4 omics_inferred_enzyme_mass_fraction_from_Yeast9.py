# Note:
# all these analysis is based on the protein abundance in the unit of mmol/DCW

# import self function
from src.protein_process import *

# data preprocess
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")

"""
# test
# omics_combine_auto.columns = omics_combine_auto.columns.str.replace('all_gene', 'gene')
# change the unit from mmol/gDCW into g/gDCW
# some general datasets
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

omics_combine_input = omics_combine_auto
omics_combine_input["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], omics_combine_input["gene"])
omics_combine_input = omics_combine_input[~omics_combine_input["MW_Kda"].isna()]

all_colum = omics_combine_input.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']

omics_combine_input2 = omics_combine_input[all_colum2]

# change mmol/gDW into g/gDW:
for x in all_colum2:
    omics_combine_input2[x] = omics_combine_input[x]*omics_combine_input["MW_Kda"]

omics_combine_input2['gene'] = omics_combine_input['gene']


# calcualte the mass fraction of enzyme under different conditions
gene_list_in_ETFL = pd.read_excel("data/proteomics/gene_list_in_ETFL.xlsx")
omics_combine_input2_for_ETFL = omics_combine_input2[omics_combine_input2['gene'].isin(gene_list_in_ETFL['geneID'])]
total_enzyme = omics_combine_input2_for_ETFL.sum(numeric_only=True, axis=0)
total_proteome = omics_combine_input2.sum(numeric_only=True, axis=0)
total_enzyme.to_excel("data/proteomics/enzyme_mass_fraction_under_each_condition.xlsx")
total_proteome.to_excel("data/proteomics/proteome_mass_faction_under_each_condition.xlsx")
"""

gene_list_in_ETFL = pd.read_excel("data/proteomics/yeast-GEM.xlsx")
# calculate mass ratio of enzyme per total protein from each organelle
def enzyme_per_protein_from_organelle(protein_abundance, enzyme_list=gene_list_in_ETFL, compartment_type="organelle"):
    """
    This function is used to calculate the mass fraction of enzyme per total protein for each organelle
    :param protein_abundance (the unit is g/gDW):
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(type="all")  # based on the automatic way
        compartment = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
        all_compartment = list(compartment.keys())

    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    # creat a dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            # test
            # y = "cytosol"
            pro_abundance = protein_abundance[['gene', col0]]
            pro_abundance.columns = ['gene', 'g/gDW']
            genes_select = compartment[y]
            # get the sum
            pro_abundance.fillna(0, axis=1, inplace=True)
            pro_select = pro_abundance[pro_abundance['gene'].isin(genes_select)]
            pro_select_enzyme = pro_select[pro_select['gene'].isin(enzyme_list["geneID"])]
            sum_all = sum(pro_abundance['g/gDW'])
            sum_select = sum(pro_select['g/gDW'])
            sum_enzyme = sum(pro_select_enzyme['g/gDW'])
            if sum_select > 0:
                ratio = sum_enzyme/sum_select
                value1.append(ratio)
            else:
                value1.append(None)
        result1[col0] = value1
    return result1

# test the above code
out1 = enzyme_per_protein_from_organelle(protein_abundance=omics_combine_auto, enzyme_list=gene_list_in_ETFL, compartment_type="organelle")
out1.to_excel("data/proteomics/enzyme_per_protein_across_compartment.xlsx")




# calculate mass ratio of enzyme per total protein within cell
def enzyme_per_total_protein_across_organelle(protein_abundance, enzyme_list=gene_list_in_ETFL, compartment_type="organelle"):
    """
    This function is used to calculate the mass fraction of enzyme per total protein for each organelle
    :param protein_abundance (the unit is g/gDW):
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        if compartment_type == "organelle":
            # compartment info
            compartment = getCompartmentGeneList(type="all")  # based on the automatic way
            compartment = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
            all_compartment = list(compartment.keys())

    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    # creat a dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            # test
            # y = "cytosol"
            pro_abundance = protein_abundance[['gene', col0]]
            pro_abundance.columns = ['gene', 'g/gDW']

            genes_select = compartment[y]
            # get the sum
            pro_abundance.fillna(0, axis=1, inplace=True)
            pro_select = pro_abundance[pro_abundance['gene'].isin(genes_select)]
            pro_select_enzyme = pro_select[pro_select['gene'].isin(enzyme_list["geneID"])]
            sum_all = sum(pro_abundance['g/gDW'])  # sum of all proteins
            sum_select = sum(pro_select['g/gDW'])  # sum of proteins from specific compartment
            sum_enzyme = sum(pro_select_enzyme['g/gDW']) # sum of proteins belong to enzymes from specific comparment
            ratio = sum_enzyme/sum_all
            value1.append(ratio)
        result1[col0] = value1
    return result1
# test the above code
out2 = enzyme_per_total_protein_across_organelle(protein_abundance=omics_combine_auto, enzyme_list=gene_list_in_ETFL, compartment_type="organelle")
out2.to_excel("data/proteomics/enzyme_fraction_based_on_total_protein_across_compartment.xlsx")
