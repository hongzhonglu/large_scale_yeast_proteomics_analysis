import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# Part 2 Collect all the data in the molecular/cell
# Generally, there are three sources, SGD, cell system and another paper. The SGD data use the median value from cell system.
# input data from SGD
# pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')

# input data from cell system, 2018
pro_abundance = pd.read_excel("data/ProteomicsData_cell_systems_2018/mmc5_normal_conditions.xlsx")
col0 = list(pro_abundance.columns)
col0 = [x for x in col0 if x != "Coefficient of Variation"]
pro_abundance = pro_abundance[col0]
col0 = col0[0:4] + [x+"_ref" for x in col0[4:]]
pro_abundance.columns = col0

pro_abundance1 = pd.read_excel("data/ProteomicsData_cell_systems_2018/mmc9_normal_plus_stress_conditions.xlsx")


# input data from cell reports 2017
# note this data is obtained under exponential growth phases
pro_abundance2 = pd.read_excel("data/proteomics/protein_copy_cell_report_2017.xlsx")
# filter out one sample with very few total protein copy number
pro_abundance2 =pro_abundance2[[x for x in pro_abundance2.columns if x !="Chong et al. 2015 - Copy "]]



# combine data from different source?
protein_copy = pd.merge(left=pro_abundance, right=pro_abundance1, left_on=['Systematic Name'], right_on=['Systematic Name'], how="outer")
protein_copy1 = pd.merge(left=protein_copy, right=pro_abundance2, left_on=['Systematic Name'], right_on=['gene'], how="outer")
col2 = list(protein_copy1.columns)
col2 = [x for x in col2 if "_x" not in x]
col2 = [x for x in col2 if "_y" not in x]
col2 = [x for x in col2 if x !="gene"]
protein_copy1 = protein_copy1[col2]
protein_copy1 = protein_copy1.rename(columns={'Systematic Name': 'gene'})

protein_copy1.to_excel("data/proteomics/protein_copy_combine.xlsx",index=False)






# change the data as mass fraction for each protein per gram of total protein mass
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")


protein_copy1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], protein_copy1["gene"])

protein_copy1 = protein_copy1[~protein_copy1["MW_Kda"].isna()]

all_colum = protein_copy1.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']

omics_combine_input2 = protein_copy1[all_colum2]


# change copy/cell into g/g total protein:
for x in all_colum2:
    omics_combine_input2[x] = protein_copy1[x]*protein_copy1["MW_Kda"]
omics_combine_input_mass_fraction = omics_combine_input2

for x in all_colum2:
    omics_combine_input_mass_fraction[x] = omics_combine_input2[x]/omics_combine_input2[x].sum()

omics_combine_input_mass_fraction["gene"] = protein_copy1["gene"]
new_column = ["gene"] + all_colum2
mass_fraction_cell_system_2018 = omics_combine_input_mass_fraction[new_column]



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
out = ProMassRatio_Organelle(protein_abundance=mass_fraction_cell_system_2018, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment_cell_system_2018.xlsx")



# how to further calculation the protein volume ratio and protein area ratio of main organelle
# change the unit from g/g into mol/g?
mass_fraction_cell_system_2018["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], mass_fraction_cell_system_2018["gene"])

all_colum = mass_fraction_cell_system_2018.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']
protein_in_mol = mass_fraction_cell_system_2018[all_colum2]
# change copy/cell into g/g total protein:
for x in all_colum2:
    protein_in_mol[x] = 1000 * protein_in_mol[x] / mass_fraction_cell_system_2018["MW_Kda"]

protein_in_mol["gene"] = mass_fraction_cell_system_2018["gene"]
new_column = ["gene"] + all_colum2

protein_in_mol = protein_in_mol[new_column]

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
            else:
                x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
                value1.append(x/total_volume)
                #value2.append(S)
        result1[col0] = value1
    return result1

s2 =Pro_3D_Volume_Ratio_Cal(protein_copy=protein_in_mol, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_cell_system_2018.xlsx")



# calculate the membrane ratio
# all datasets
s2 =ProMembraneCal(protein_in_mol) # from part 3.9
s2.to_excel("data/proteomics/membrane_size_ratio_across_compartment_cell_system_2018.xlsx")

