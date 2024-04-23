# this script is to process compartment datasets
# 2021-11-16

# here the compartment annotation is mainly from SGD and MitoMiner
import pandas as pd
from src.protein_process import *
from matplotlib import pyplot as plt

# test the function
compartment_dict_all0 = getCompartmentGeneList(type="all")
compartment_dict20 = getCompartmentGeneList(type="manual") # Remove some compartmental annotation only with computational evidence (keep experimental evidence)


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
compartment_compare.columns = ["compartment", "all_annotation", "annotation_manual_evidance"]
compartment_compare.to_excel("data/compare_compartment_anotation_with_all_and_manual_evidance.xlsx")


# add one plot
df2_order = df2.sort_values(by=['gene_number'], ascending=False)
x_value = df2_order["compartment"].tolist()[0:10]
y_value = df2_order["gene_number"].tolist()[0:10]
plt.bar(x_value, y_value, width=0.7, bottom=50, align='edge')
plt.title('')
plt.xlabel('Compartment', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(rotation=90)
plt.show()

