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
    return coef[0], coef[1]

structure_quality_all = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
pro_size00 = pro_size[pro_size["locus"].isin(structure_quality_all["gene"])]

pro_size00["MW"] = pro_size00["MW"]/1000 # change the unit as kda
a, b = linearFit(df=pro_size00, x_name="MW", y_name="Total_Volume")

#pro_size00["calculated"] = pro_size00["MW"]*1.06e-03 - 1.10
pro_size00["calculated"] = pro_size00["MW"]*a - b

pro_size00["relative_change"] = (pro_size00["Total_Volume"] - pro_size00["calculated"])/pro_size00["calculated"]

# filter by protein length to remove too short proteins
# pro_size00 = pro_size00[pro_size00["pro_length"] >=200]
pro_size00 = pro_size00.sort_values(by=['relative_change'], ascending=True)


pro_size01 = pro_size00.iloc[0:200,:]
gene01= ",".join(pro_size01["locus"].to_list())
print(gene01)




# import Transcription factor
from scipy.stats import ttest_ind
TRN_sce = pd.read_excel("data/transcriptional_network/TRN_sce.xlsx")
pro_size00["TF"] = None
pro_size00["TF"][pro_size00["locus"].isin(TRN_sce["TF"])] = "Yes"
pro_size00["TF"][~pro_size00["locus"].isin(TRN_sce["TF"])] = "No"
pro_size00["volume_per_kda"] = pro_size00["Total_Volume"]/pro_size00["MW"]
pro_size00["id"] = singleMapping(structure_quality_all["id"],structure_quality_all["gene"],pro_size00["locus"])
pro_size00.to_excel("data/sce_protein_with_TF_classification.xlsx")
pro_g1 = pro_size00[pro_size00["TF"]=="Yes"]
pro_g2 = pro_size00[pro_size00["TF"]=="No"]
pro_g2 = pro_g2[pro_g2["pro_length"] >= min(pro_g1["pro_length"])]
pro_g2 = pro_g2[pro_g2["pro_length"] <= max(pro_g1["pro_length"])]
# combine two pandas
pro_c = pd.concat([pro_g1, pro_g2], axis=0)
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])



pro_size00 = pro_size00.sort_values(by=['volume_per_kda'], ascending=True)
sns.displot(pro_size00, x="volume_per_kda", stat="density", common_norm=False)
plt.xlim(0.75,1.25)
pro_size01 = pro_size00.iloc[0:200,:]
gene01= ",".join(pro_size01["locus"].to_list())
print(gene01)




# compare new and conserved gene
gene_type = pd.read_excel("data/core_and_essential_gene_sce/panGene_for manual check.xlsx")
pro_size00["gene_type"] = singleMapping(gene_type["gene_type"],gene_type["gene_simple"],pro_size00["locus"])
pro_size00["volume_per_kda"] = pro_size00["Total_Volume"]/pro_size00["MW"]
pro_size00 = pro_size00[~pro_size00["gene_type"].isna()]
pro_size00["id"] = singleMapping(structure_quality_all["id"],structure_quality_all["gene"],pro_size00["locus"])
pro_size00.to_excel("data/sce_protein_with_core_variable_type.xlsx")
sns.catplot(x="gene_type", y="volume_per_kda", order=["core_gene", "Variable"], kind="box", data=pro_size00)
pro_g1 = pro_size00[pro_size00["gene_type"]=="core_gene"]
pro_g2 = pro_size00[pro_size00["gene_type"]=="Variable"]
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])
pro_g1['volume_per_kda'].describe()
pro_g2['volume_per_kda'].describe()




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


# based on quality, estimate which enzyme from Yeast8 need re-modelling
data_merge = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
sce_pathway["id"] = singleMapping(data_merge["id"],data_merge["gene"],sce_pathway["gene"])
sce_pathway.to_excel("data/sce_kegg_pathway.xlsx")




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

