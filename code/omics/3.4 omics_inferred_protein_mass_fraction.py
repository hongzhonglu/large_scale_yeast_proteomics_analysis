# Note:
# all these analysis is based on the protein abundance in the unit of mmol/DCW


# import self function
from src.model_process import *
from src.protein_process import *

# data preprocess
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx") # the unit the g/gDW???? should be wrong

# test
omics_combine_auto.columns = omics_combine_auto.columns.str.replace('all_gene', 'gene')
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
for x in all_colum2:
    omics_combine_input2[x] = omics_combine_input[x]*omics_combine_input["MW_Kda"]

omics_combine_input2['gene'] = omics_combine_input['gene']


# calculate mass ratio of protein from each organelle per total protein mass
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
out = ProMassRatio_Organelle(protein_abundance=omics_combine_input2, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment.xlsx")

# calcualte the mass fraction of enzyme under different conditions
gene_list_in_ETFL = pd.read_excel("data/proteomics/gene_list_in_ETFL.xlsx")
omics_combine_input2_for_ETFL = omics_combine_input2[omics_combine_input2['gene'].isin(gene_list_in_ETFL['geneID'])]
total_enzyme = omics_combine_input2_for_ETFL.sum(numeric_only=True, axis=0)
total_proteome = omics_combine_input2.sum(numeric_only=True, axis=0)
total_enzyme.to_excel("data/proteomics/enzyme_mass_fraction_under_each_condition.xlsx")
total_proteome.to_excel("data/proteomics/proteome_mass_faction_under_each_condition.xlsx")
