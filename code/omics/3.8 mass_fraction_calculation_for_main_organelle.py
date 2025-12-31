# note:
# this part can be rewritten as a function

import os
# import self function
from src.protein_process import *

compartment_all = getCompartmentGeneList(type="all")
compartment_list = list(compartment_all.keys())
compartment_out = ','.join(compartment_list)

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
    "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
    "plasma membrane", "P-body", "cytoplasmic stress granule", "ribosome","cytosol",
    "cytoskeleton","extracellular region"]


compartment_main = {k: v for k, v in compartment_all.items() if k in organelles_main}

compartment_main_c = gene_location_curation_sce(compartment_main)
compartment_main_c = {k: v for k, v in compartment_main_c.items() if k in organelles_main}



from collections import defaultdict
def count_keys_per_value_optimized(input_dict):
    # 自动初始化值为 list 类型
    result = defaultdict(list)

    for key, value_list in input_dict.items():
        for item in value_list:
            result[item].append(key)

    return dict(result)  # 转回普通字典返回

output = count_keys_per_value_optimized(compartment_main_c)
print(output)









def ProMassRatio_Organelle(protein_abundance, compartment=compartment_all):
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




import pandas as pd

Sample_ID_select = list(mass_fraction_final.columns)
Sample_ID_select = [x for x in Sample_ID_select if x !="gene"]
for col0 in Sample_ID_select:
    print(col0)
    value1 = []
    pro_abundance = mass_fraction_final[['gene', col0]]
    pro_abundance.columns = ['gene', 'g/gDW']
    pro_abundance = pro_abundance.dropna()

    gene = pro_abundance['gene'].tolist()
    gene_with_compartment = list(output.keys())
    gene_refine = list(set(gene) & set(gene_with_compartment))

    pro_abundance00 = pro_abundance[pro_abundance['gene'].isin(gene_refine)]
    gene_refine00 = pro_abundance00['gene'].tolist()

    pro_local = []
    for x in gene_refine00:
        print(list(output[x]))
        pro_local.append(list(output[x]))

    # 示例数据（实际替换为您的数据集，如从SGD或Ho et al.下载）
    data = {
        'protein': pro_abundance00['gene'].tolist(),
        'abundance': pro_abundance00['g/gDW'].tolist(),  # 每细胞拷贝数
        'localizations': pro_local
    }
    df = pd.DataFrame(data)

    # 您的细胞器列表
    organelles_main = [
        "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
        "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
        "plasma membrane", "P-body", "cytoplasmic stress granule", "ribosome", "cytosol",
        "cytoskeleton", "extracellular region"
    ]

    # Fractional分配计算
    allocation = {org: 0.0 for org in organelles_main}  # 初始化
    total_abundance = df['abundance'].sum()  # 总资源（守恒）

    for _, row in df.iterrows():
        locs = [org for org in row['localizations'] if org in organelles_main]  # 只计您的列表中
        if locs:  # 避免未注解蛋白
            fraction = row['abundance'] / len(locs)  # 均匀分配
            for loc in locs:
                allocation[loc] += fraction

    # 计算比例
    proportions = {org: (allocation[org] / total_abundance) * 100 for org in organelles_main}

    # 输出（排序显示）
    proportions_sorted = dict(sorted(proportions.items(), key=lambda x: x[1], reverse=True))
    for org, pct in proportions_sorted.items():
        print(f"{org}: {pct:.2f}%")









