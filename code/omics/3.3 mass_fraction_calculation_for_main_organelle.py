# note:
# this part can be rewritten as a function
import pandas as pd
import os
# import self function
from src.protein_process import *

##
mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")

compartment_all0 = getCompartmentGeneList(type="all")

compartment_all00 = gene_location_curation_sce(organelle0=compartment_all0)

## 重新整合，包括实验、手动校正和所有
# 这个23个细胞器或者子细胞器完全基于有实验证据的，且未经手工查询校正
Exp_list = pd.read_excel('data/compare_compartment_anotation_from_different_version_check.xlsx', sheet_name='g2-exp')
compartment_dict20 = getCompartmentGeneList(type="manual") # Remove some compartmental annotation only with computational evidence (keep experimental evidence)
filtered = {k: v for k, v in compartment_dict20.items() if k in Exp_list['compartment'].tolist()}

# 其他103个， 21个细胞器是经过手动校正，剩余的是基于全部证据整合(其中就包括计算推测与实验)
filtered2 = {k: v for k, v in compartment_all00.items() if k not in Exp_list['compartment'].tolist()}
# 去重
duplicated = ['mitochondrial membrane','vacuolar membrane','vacuole']
filtered20 = {k: v for k, v in filtered2.items() if k not in duplicated}
# 合并
compartment_all = {**filtered, **filtered20}



compartment_list = list(compartment_all.keys())
compartment_out = ','.join(compartment_list)

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
    "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
    "plasma membrane", "P-body", "spindle pole body", "ribosome", "cytosol", "extracellular region"]
# "cytosol","extracellular region" 属于区室，而非细胞器


compartment1 = {k: v for k, v in compartment_all.items() if k in organelles_main}
# remove duplicates
compartment1 = {key: list(set(value)) for key, value in compartment1.items()}


mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")


out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment1)
out00.to_excel("data/proteomics/main_organelle_fraction_test.xlsx")


out11 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment_all)
out11.to_excel("data/proteomics/all_organelle_fraction_test.xlsx")


# 逐行计算统计量
df = out00
df_stats = pd.DataFrame()
# 计算均值
df_stats['mean'] = df.mean(axis=1)
# 计算标准差
df_stats['std'] = df.std(axis=1)
# 计算分位数
df_stats['q25'] = df.quantile(0.25, axis=1)
df_stats['q50'] = df.quantile(0.5, axis=1)  # 中位数
df_stats['q75'] = df.quantile(0.75, axis=1)
# 计算最小值和最大值
df_stats['min'] = df.min(axis=1)
df_stats['max'] = df.max(axis=1)
# 计算范围（极差）
df_stats['range'] = df.max(axis=1) - df.min(axis=1)
df_stats.index = df['compartment'].tolist()