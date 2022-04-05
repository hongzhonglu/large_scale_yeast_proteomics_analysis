# gene level analysis

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")

# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]




# further input the protein length and molecular weight
pro_info = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")

# further input the subsystem information of each gene
sub_pathway = pd.read_csv("data/subsystem/pathway_list_kegg.txt", header=None, sep="\t")
sub_pathway.columns = ["pathway","name"]
sce_pathway = pd.read_csv("data/subsystem/sce_pathway.txt", header=None, sep="\t")
sce_pathway.columns = ["pathway","gene"]
sce_pathway["name"] = singleMapping(sub_pathway["name"],sub_pathway["pathway"],sce_pathway["pathway"])
sce_pathway["gene"] = sce_pathway["gene"].str.replace("sce:", "")

sce_pathway["MW"] = singleMapping(pro_info["proteins_molecular_weight"], pro_info["locus"],sce_pathway["gene"])
sce_pathway["pro_length"] = singleMapping(pro_info["protein_length"],pro_info["locus"],sce_pathway["gene"])
sce_pathway["pro_volume"] = singleMapping(pro_size["Total_Volume"],pro_size["locus"],sce_pathway["gene"])

# two group
pathway_select = ['Glycolysis / Gluconeogenesis', 'Citrate cycle (TCA cycle)','Pentose phosphate pathway']
sce_pathway["group"] = None
sce_pathway["group"][sce_pathway["name"].isin(pathway_select)] = "core"
sce_pathway["group"][~sce_pathway["name"].isin(pathway_select)] = "other"

# Remove duplicates
sce_pathway1 = sce_pathway[sce_pathway["group"]=="core"]
sce_pathway1 = sce_pathway1.drop_duplicates('gene', keep='last')
sce_pathway2 = sce_pathway[sce_pathway["group"]=="other"]
sce_pathway2 = sce_pathway2.drop_duplicates('gene', keep='last')
sce_pathway_update = pd.concat([sce_pathway1, sce_pathway2])
# protein volume over 100nm^3
sce_pathway0 = sce_pathway[sce_pathway["pro_volume"] >100]


sns.displot(sce_pathway_update, x="pro_volume", hue="group", stat="density", common_norm=False)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)




# complex  or non complex??
complex_inf = pd.read_excel("data/complex_info.xlsx")
complex_inf["subunit"] = complex_inf["subunit"].str.replace("-MONOMER","")
pro_size["complex"] = None
pro_size["complex"][pro_size["locus"].isin(complex_inf["subunit"])] = "is_complex"
pro_size["complex"][~pro_size["locus"].isin(complex_inf["subunit"])] = "not_complex"
sns.displot(pro_size, x="Total_Volume", hue="complex", stat="density", common_norm=False)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# for all the proteins
sns.displot(pro_size, x="Total_Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["Total_Volume"]]
df1 = df.describe()

sns.displot(pro_size, x="section_area_new")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["section_area_new"]]
df2 = df.describe()

