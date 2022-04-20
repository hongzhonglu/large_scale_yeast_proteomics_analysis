# this script is to process compartment datasets
# 2021-11-16

# here the compartment annotation is mainly from SGD and MitoMiner
import pandas as pd
from src.protein_process import *


# test the function
compartment_dict_all0 = getCompartmentGeneList(filter="No")
compartment_dict20 = getCompartmentGeneList(filter="Yes")


# compare the difference
key0 = []
len0 = []

for key in compartment_dict_all0.keys():
    key0.append(key)
    len0.append(len(compartment_dict_all0[key]))
df1 = pd.DataFrame({"compartment":key0, "gene_number": len0})

key0 = []
len0 = []
for key in compartment_dict20.keys():
    key0.append(key)
    len0.append(len(compartment_dict20[key]))
df2 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# combine the dataframe
compartment_compare = pd.merge(left=df1, right=df2, left_on=['compartment'], right_on=['compartment'], how='left')
compartment_compare.columns = ["compartment", "all_annotation", "annotation_filter"]
compartment_compare.to_excel("data/compare_compartment_anotation_with_and_without_filter.xlsx")