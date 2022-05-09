# gene level analysis

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the pro structure size data
pro_size = pd.read_excel("data/other_species/ecoli_structure_size_part1.xlsx")

pro_size = pro_size[['Protein','Total_Volume']]
pro_size["Protein"] = pro_size["Protein"].str.replace("-F1-model_v2","").str.replace("AF-","")

# for all the proteins
sns.displot(pro_size, x="Total_Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["Total_Volume"]]
df1 = df.describe()

# further input the protein length and molecular weight
pro_info = pd.read_csv("data/other_species/ecoli.txt", sep="\t")
pro_size["MW"] = singleMapping(pro_info["Mass"], pro_info["Entry"],pro_size["Protein"])
pro_size["pro_length"] = singleMapping(pro_info["Length"], pro_info["Entry"],pro_size["Protein"])

pro_size["MW"] = pro_size["MW"].str.replace(",", "")
pro_size.MW = pd.to_numeric(pro_size.MW, errors='coerce')


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

a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
# filter by protein length to remove too short proteins

pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]
pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)









# spo
# input the pro structure size data
pro_size = pd.read_excel("data/other_species/spo_structure_size_part1.xlsx")

pro_size = pro_size[['Protein','Total_Volume']]
pro_size["Protein"] = pro_size["Protein"].str.replace("-F1-model_v2","").str.replace("AF-","")

# for all the proteins
sns.displot(pro_size, x="Total_Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
df = pro_size[["Total_Volume"]]
df1 = df.describe()

# further input the protein length and molecular weight
pro_info = pd.read_csv("data/other_species/spo.txt", sep="\t")
pro_size["MW"] = singleMapping(pro_info["Mass"], pro_info["Entry"],pro_size["Protein"])
pro_size["pro_length"] = singleMapping(pro_info["Length"], pro_info["Entry"],pro_size["Protein"])

pro_size["MW"] = pro_size["MW"].str.replace(",", "")
pro_size.MW = pd.to_numeric(pro_size.MW, errors='coerce')


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

a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
# filter by protein length to remove too short proteins

pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]
pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)

sns.displot(pro_size, x="volume_per_kda", stat="density", common_norm=False)


# import Transcription factor
from scipy.stats import ttest_ind
pro_size00 = pro_size.copy()
TRN_spo = pd.read_csv("data/other_species/spo_TF_list.tsv", sep="\t")
pro_size00["TF"] = None
pro_size00["TF"][pro_size00["Protein"].isin(TRN_spo["UniProt ID"])] = "Yes"
pro_size00["TF"][~pro_size00["Protein"].isin(TRN_spo["UniProt ID"])] = "No"
pro_g1 = pro_size00[pro_size00["TF"]=="Yes"]
pro_g2 = pro_size00[pro_size00["TF"]=="No"]
pro_g2 = pro_g2[pro_g2["pro_length"] >= min(pro_g1["pro_length"])]
pro_g2 = pro_g2[pro_g2["pro_length"] <= max(pro_g1["pro_length"])]
# combine two pandas
pro_c = pd.concat([pro_g1, pro_g2], axis=0)
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])




