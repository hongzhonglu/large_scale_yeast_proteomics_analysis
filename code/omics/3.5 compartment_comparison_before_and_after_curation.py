# note:
# this part can be rewritten as a function
import pandas as pd
import os
# import self function
from src.protein_process import *
# reanalyze the data set
# original compartment
compartment_all0 = getCompartmentGeneList(type="all")  # based on the automatic way
# remove duplicates
compartment = {key: list(set(value)) for key, value in compartment_all0.items()}
key0 = []
len0 = []
for key in compartment.keys():
    key0.append(key)
    len0.append(len(compartment[key]))
df1 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# after curation
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

key0 = []
len0 = []
for key in compartment_all.keys():
    key0.append(key)
    len0.append(len(compartment_all[key]))
df2 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# combine the dataframe
compartment_compare = pd.merge(left=df1, right=df2, left_on=['compartment'], right_on=['compartment'], how='outer')
compartment_compare.columns = ["compartment", "annotation_combine", "annotation_curation"]
compartment_compare.to_excel("data/compare_compartment_annotation_with_and_without_manual_curation.xlsx")





