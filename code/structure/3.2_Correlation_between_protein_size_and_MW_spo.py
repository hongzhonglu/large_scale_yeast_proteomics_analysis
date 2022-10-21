# gene level analysis

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from scipy.stats import ttest_ind
# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *


def calculateQualityScore(pdb_dir="/Users/xluhon/Documents/alphafold_pdb/"):
    import os
    import pandas as pd
    from biopandas.pdb import PandasPdb
    pdbfile = pdb_dir
    all_pdb = os.listdir(pdbfile)
    ppdb = PandasPdb()
    PITT_score = []
    for i in all_pdb:
        print(i)
        if '.pdb' in i:
            pdb_in = pdbfile + i
            ppdb.read_pdb(pdb_in)
            ss = ppdb.df['ATOM']
            mainchain = ss[(ss['atom_name'] == 'CA')]
            bfact_mc_avg = mainchain['b_factor'].mean()
            PITT_score.append(bfact_mc_avg)
    data_merge = pd.DataFrame({"id": all_pdb, "score": PITT_score})
    data_merge["id_update"] = data_merge["id"].str.replace("AF-", "").str.replace("-F1-model_v1.pdb", "")
    return data_merge



# spo
# lack of quality score
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
quality_spo = calculateQualityScore(pdb_dir='/Users/xluhon/Documents/alphafold_pdb_SCHPO_v2/')

quality_spo['id_update'] = quality_spo['id_update'].str.replace('-F1-model_v2.pdb','')

pro_size["score"] = singleMapping(quality_spo['score'], quality_spo['id_update'], pro_size["Protein"])



# quality analysis as a whole
score_list = quality_spo["score"].tolist()
score_high =[x for x in score_list if x >= 75] # 60.66% proteins are of high-quality
print(len(score_high)/len(score_list))


sns.displot(quality_spo, bins=14, x="score",alpha=.4, height=3, aspect=1.2)
plt.xlabel('pLDDT average score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/structure_quality_spo.pdf', bbox_inches='tight')








a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]
sns.displot(pro_size, x="volume_per_kda", stat="density", common_norm=False)

# enrichment analysis
pro_size = pro_size.dropna()
pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)

pro_size02 = pro_size.iloc[4912:5112,:]
gene02= ",".join(pro_size02["Protein"].to_list())
print(gene02)




# one interesting idea is used sce formula to calculate the protein volume in other species
pro_size["calculated_volume"] = 1.06019171e-03*pro_size["MW"] - 1.10587455
# compare the predicted and calculated for e.coli
x0 = "MW"
y0 = "Total_Volume"
y1 = "calculated_volume"
plt.figure(figsize=(3, 3.6))
sns.lineplot(x=x0, y=y1, data=pro_size)
sns.scatterplot(x=x0, y=y0, data=pro_size)
plt.xlabel(x0, fontsize=15)
plt.ylabel(y1, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/fitted_structure_volume_spo.pdf', bbox_inches='tight')







# import Transcription factor
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
sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=pro_c)



pro_c.to_excel('data/other_species/spo_structure_info.xlsx')



# re-do the above analysis using the calibrated datasets
spo_calibrated = pd.read_excel('data/other_species/spo_structure_info2.xlsx')
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=spo_calibrated)

sns.catplot(x="TF", y="volume_per_kda2", order=["No", "Yes"], kind="box", data=spo_calibrated)
plt.xlabel("TF",fontsize=15)
plt.ylabel("Volume_per_kda",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=spo_calibrated)
pro_g1 = spo_calibrated[spo_calibrated["TF"]=="Yes"]
pro_g2 = spo_calibrated[spo_calibrated["TF"]=="No"]
ttest_ind(pro_g1['volume_per_kda2'], pro_g2['volume_per_kda2'])


spo_calibrated2 = spo_calibrated[spo_calibrated['score'] >=75]
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=spo_calibrated2)
sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=spo_calibrated2)
pro_g1 = spo_calibrated2[spo_calibrated2["TF"]=="Yes"]
pro_g2 = spo_calibrated2[spo_calibrated2["TF"]=="No"]
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])
