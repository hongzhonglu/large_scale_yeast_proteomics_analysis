# note:
# this part can be rewritten as a function

import os
# import self function
import pandas as pd
from src.protein_process import *

mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")


compartment_all0 = getCompartmentGeneList(type="all")
compartment_all = gene_location_curation_sce(organelle0=compartment_all0)

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


# quality check
plasma = compartment_main_c['plasma membrane']
mass_fraction_final_plasma = mass_fraction_final[mass_fraction_final['gene'].isin(plasma)]
column_sums1 = mass_fraction_final_plasma.sum()
column_sums2 = mass_fraction_final.sum()



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











# new way
Sample_ID_select = list(mass_fraction_final.columns)
Sample_ID_select = [x for x in Sample_ID_select if x !="gene"]
output0 = {}
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
    output0[col0] = proportions_sorted



output_df = pd.DataFrame(output0)
output_df = output_df / 100
output_df.to_excel("data/sce_compartment_curation/mass_fraction_new_way.xlsx")




