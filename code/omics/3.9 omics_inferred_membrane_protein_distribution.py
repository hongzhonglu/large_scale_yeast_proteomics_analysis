# Note:
# which protein on the membrane
# 1) all the transporter protein
# 2) all the membrane-related protein


# import self function
from src.model_process import *
from src.protein_process import *


# as the first step: define the membrane or transporter protein
protein_transporter = open("/Users/xluhon/Documents/GitHub/large_scale_yeast_proteomics_analysis/data/tcdb.txt").readlines()
protein_transporter = [x for x in protein_transporter if ">" in x]
protein_transporter = [x for x in protein_transporter if "S288c" in x]
protein_ID = []
for xx in protein_transporter:
    ss0 = xx.split("|")[2]
    protein_ID.append(ss0)
# get the transporter gene id in sce
uniprotGeneID_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
transporter_tf = uniprotGeneID_mapping[uniprotGeneID_mapping["Entry"].isin(protein_ID)]
transporter_pro_list = transporter_tf["GeneName"].tolist()

# get the membrane annotation from SGD
membrane_pro = pd.read_excel("data/membrane_annotations.xlsx")
membrane_pro_list = list(set(membrane_pro['Systematic Name/Complex Accession'].tolist()))

# Input the datasets from paxDB
compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')
# extract compartment
compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier',
                       'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]
compartment1_membrane_filter = compartment1[compartment1["GO_Name"].str.contains("membrane")]
membrane_pro_list_database = list(set(compartment1_membrane_filter["Systematic_name"].tolist()))


# check the relation between the annotation from the above procedures
membrane_pro_final_merge = list(set(membrane_pro_list) & set(membrane_pro_list_database))


# update this function
def ProMembraneCal(protein_copy, compartment_type="organelle"):
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
    all_compartment = [x for x in all_compartment if "membrane" in x]
    result2 = pd.DataFrame({"compartment": all_compartment})
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



# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")


# only analyze the rosemary datasets under NH4 limitation?
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all_rosemary = protein_copy_all1[Sample_ID_select+["gene"]]
# curate the protein copies based on newly calculated size
df_curated = calculateCurationCoefficent()
curation_coefficent = df_curated["curation_coefficent"].to_list()
for i, sid in enumerate(Sample_ID_select):
    print(i, sid, curation_coefficent[i])
    protein_copy_all_rosemary[sid] = protein_copy_all_rosemary[sid]*curation_coefficent[i]

protein_copy_all_rosemary.to_excel("data/proteomics/all_protein_copy_rosemary.xlsx")

s2 = ProMembraneCal(protein_copy_all_rosemary)
s2.to_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v3.xlsx")


# jianye datasets
all_columns = list(protein_copy_all1.columns)
Sample_ID_select = [x for x in all_columns if "_M" in x]
protein_copy_jianye = protein_copy_all1[Sample_ID_select+["gene"]]
s2 = ProMembraneCal(protein_copy_jianye)
s2.to_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")


# all datasets
s2 =ProMembraneCal(protein_copy_all1)
s2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")
