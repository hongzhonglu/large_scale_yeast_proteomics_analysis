# note:
# this part can be rewritten as a function

import os
# import self function
from src.protein_process import *

"""
# Part 1
# Initially check how many genes could find compartment
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
#compartment = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
all_compartment = list(compartment.keys())
mitochondrion_related = ['mitochondrion',
'mitochondrial outer membrane',
'mitochondrial intermembrane space',
'mitochondrial matrix',
'mitochondrial inner membrane']
compartment_m = {x:y for x, y in compartment.items() if x in mitochondrion_related}
"""


# reanalyze the data set
gene_mitochondrion_sgd = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx") # this is based sgd annotation
mitochondrion_gene = gene_mitochondrion_sgd['gene'].tolist()

# combine
mitochondrial_inner_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_outer_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_matrix = pd.read_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v2.xlsx")

# combine
mitochondrial_intermembrane_space_annotations = pd.read_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v2.xlsx")




# input the NC paper
gene_mitochondrion_nc = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
gene_mitochondrion_nc = gene_mitochondrion_nc.iloc[0:987:]
gene_m_Outer_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Outer membrane'] == "X"]
gene_m_Inner_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Inner membrane'] == "X"]
gene_m_OI_space = gene_mitochondrion_nc[gene_mitochondrion_nc['Inter-membrane space'] == "X"]
gene_m_matrix = gene_mitochondrion_nc[gene_mitochondrion_nc['Matrix'] == "X"]



gene_m_Outer_membrane = list(set(gene_m_Outer_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_Inner_membrane = list(set(gene_m_Inner_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_OI_space = list(set(gene_m_OI_space['gene'].tolist()) & set(mitochondrion_gene))
gene_m_matrix = list(set(gene_m_matrix['gene'].tolist()) & set(mitochondrion_gene))



# further combine different data together:
gene_m_Outer_membrane_v3 = list(set(gene_m_Outer_membrane + mitochondrial_outer_membrane['gene'].tolist()))
gene_m_Inner_membrane_v3 = list(set(gene_m_Inner_membrane + mitochondrial_inner_membrane['gene'].tolist()))
gene_m_OI_space_v3 = list(set(gene_m_OI_space + mitochondrial_intermembrane_space_annotations['gene'].tolist()))
gene_m_matrix_v3 = list(set(gene_m_matrix + mitochondrial_matrix['gene'].tolist()))
gene_assigned = list(set(gene_m_Outer_membrane_v3 + gene_m_Inner_membrane_v3 + gene_m_OI_space_v3 + gene_m_matrix_v3)) # 742 proteins
gene_unassigned = list(set(mitochondrion_gene)-set(gene_assigned)) # 454 proteins



gene_m_Outer_membrane_v3 = pd.DataFrame({"gene": gene_m_Outer_membrane_v3})
gene_m_Inner_membrane_v3 = pd.DataFrame({"gene": gene_m_Inner_membrane_v3})
gene_m_OI_space_v3 = pd.DataFrame({"gene": gene_m_OI_space_v3})
gene_m_matrix_v3 = pd.DataFrame({"gene": gene_m_matrix_v3})
gene_m_unassigned_v3 = pd.DataFrame({"gene": gene_unassigned})


gene_m_Outer_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v3.xlsx")
gene_m_Inner_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v3.xlsx")
gene_m_OI_space_v3.to_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v3.xlsx")
gene_m_matrix_v3.to_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v3.xlsx")
gene_m_unassigned_v3.to_excel("data/sce_compartment_curation/mitochondrial_unassigned_manual_v3.xlsx")


# input mitochrondrial protein list from  Absolute yeast mitochondrial proteome quantification reveals trade-off between biosynthesis and energy generation during diauxic shift
# input the data from PNAS
gene_mitochondrion_pnas = pd.read_excel("data/sce_compartment_curation/pnas.1918216117.sd01 VIP Francesca.xlsx", header=1)
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
sce_gene_select =sce_gene[sce_gene['Entry'].isin(gene_mitochondrion_pnas['Uniprot Accession'])]
gene_mitochondrion_pnas_F = sce_gene_select

# input mitochrondrial protein list from  jianye xia
gene_mitochondrion_xia = pd.read_excel("data/sce_compartment_curation/41467_2022_30513_MOESM5_ESM xia.xlsx", header=0)
gene_mitochondrion_xia.columns =['gene','short_name']
gene_mitochondrion_xia.to_excel("data/sce_compartment_curation/mitochondrion_xia.xlsx")



# input mitochrondrial protein list from  Definition of a High-Confidence Mitochondrial Proteome at Quantitative Scale
gene_mitochondrion_marcel = pd.read_excel("data/sce_compartment_curation/mmc4_cell_reports_2017_marcel.xlsx", header=0)
gene_mitochondrion_marcel = gene_mitochondrion_marcel.iloc[:,0:2]
gene_mitochondrion_marcel.columns =['gene','Entry']




compartment1 = {}
compartment1['m_nc'] = gene_mitochondrion_nc['gene'].tolist()
compartment1['m_sgd'] = gene_mitochondrion_sgd['gene'].tolist()
compartment1['m_xia'] = gene_mitochondrion_xia['gene'].tolist()
compartment1['m_pnas'] = gene_mitochondrion_pnas_F['GeneName'].tolist()
compartment1['m_marcel'] = gene_mitochondrion_marcel['gene'].tolist()




compartment_all = getCompartmentGeneList(type="all")
compartment_list = list(compartment_all.keys())
compartment_out = ','.join(compartment_list)

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
    "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
    "plasma membrane", "P-body", "cytoplasmic stress granule", "ribosome","cytosol",
    "cytoskeleton","extracellular region"]


compartment_all123 = {k: v for k, v in compartment_all.items() if k=="mitochondrion"}
compartment_all_other = {k: v for k, v in compartment_all.items() if "itochondri" not in k}
compartment_all_f = {**compartment_all123, **compartment_all_other}
from collections import defaultdict
def count_keys_per_value_optimized(input_dict):
    # 自动初始化值为 list 类型
    result = defaultdict(list)

    for key, value_list in input_dict.items():
        for item in value_list:
            result[item].append(key)

    return dict(result)  # 转回普通字典返回

output = count_keys_per_value_optimized(compartment_all_f)
print(output)

output_f = {k: v for k, v in output.items() if 'mitochondrion' in v}
# check which protein exist only in mitochondrion
output_f2 = {k: v for k, v in output_f.items() if len(v) <=2}



compartment1['m_sgd_core'] = list(output_f2.keys())



def ProMassRatio_Organelle(protein_abundance, compartment=compartment1):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_abundance:
    :param compartment_type:
    :return:
    """
    # test
    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    # use some manually checked gene compartment definion
    # gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    # gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    # creat a dataframe to save the result

    all_compartment = list(compartment.keys())
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
mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")
out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final)
out00.to_excel("data/sce_compartment_curation/mitochondrion_fraction_under_different_input.xlsx")


