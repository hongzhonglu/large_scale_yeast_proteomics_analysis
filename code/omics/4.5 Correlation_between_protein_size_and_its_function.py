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


# further input the protein length and molecular weight
pro_info = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
pro_size["MW"] = singleMapping(pro_info["proteins_molecular_weight"], pro_info["locus"],pro_size["locus"])
pro_size["pro_length"] = singleMapping(pro_info["protein_length"], pro_info["locus"],pro_size["locus"])
plt.figure()
plt.scatter(pro_size["MW"], pro_size["Total_Volume"], marker='.')
plt.figure()
plt.scatter(pro_size["pro_length"], pro_size["Total_Volume"], marker='.')


def linearFit(df, x_name, y_name):
    """
    Linear fit for two colums in a dataframe
    :param df:
    :param x_name:
    :param y_name:
    :return:
    """
    from sklearn.metrics import r2_score
    df = df[[x_name, y_name]]
    df = df.dropna()
    x = df[x_name]
    y = df[y_name]
    x_name = x_name.split("(")[0]
    y_name = y_name.split("(")[0]
    x_name = x_name.replace("/", "_per_")
    y_name = y_name.replace("/", "_per_")
    coef = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coef)
    predict = np.poly1d(coef)
    R2 = r2_score(y, predict(x))
    print(R2)
    print(coef)
    R2 = "{:.3f}".format(R2)
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.figure()
    plt.plot(x, y, 'yo', x, poly1d_fn(x), '--k')  # '--k'=black dashed line, 'yo' = yellow circle marker
    plt.xlabel(x_name, fontsize=18)
    plt.ylabel(y_name, fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    x_max = max(x)
    y_max = max(y)
    plt.text(x_max/3, 2*y_max/3, "R2=" + str(R2), fontsize=18)
    plt.show()
structure_quality_all = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
pro_size00 = pro_size[pro_size["locus"].isin(structure_quality_all["gene"])]
linearFit(df=pro_size00, x_name="MW", y_name="Total_Volume")

pro_size00["calculated"] = pro_size00["MW"]*1.06e-03-1.10
pro_size00["relative_change"] = (pro_size00["Total_Volume"] - pro_size00["calculated"])/pro_size00["calculated"]
pro_size00 = pro_size00.sort_values(by=['relative_change'], ascending=True)

pro_size01 = pro_size00.iloc[0:150,:]
gene01= ",".join(pro_size01["locus"].to_list())



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

