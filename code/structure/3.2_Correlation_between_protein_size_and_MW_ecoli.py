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
from scipy.stats import gaussian_kde
from scipy import stats

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





# input the pro structure size data
# lack of quality score
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
a, b = linearFit(df=pro_size, x_name="MW", y_name="Total_Volume")

pro_size["calculated"] = pro_size["MW"]*a-b
pro_size["relative_change"] = (pro_size["Total_Volume"] - pro_size["calculated"])/pro_size["calculated"]
pro_size["volume_per_kda"] = 1000*pro_size["Total_Volume"]/pro_size["MW"]


# calculate the quality score
quality_ecoli = calculateQualityScore(pdb_dir='/Users/xluhon/Documents/alphafold_pdb_ECOLI_v2/')

quality_ecoli['id_update'] = quality_ecoli['id_update'].str.replace('-F1-model_v2.pdb','')

pro_size["score"] = singleMapping(quality_ecoli['score'], quality_ecoli['id_update'], pro_size["Protein"])


# quality analysis as a whole
score_list = quality_ecoli["score"].tolist()
score_high =[x for x in score_list if x >= 75]
print(len(score_high)/len(score_list))


sns.displot(quality_ecoli, bins=20, x="score", alpha=.2, height=3, aspect=1.2)
plt.xlabel('pLDDT average score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.axvline(x=75)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/structure_quality_ecoli.pdf', bbox_inches='tight')







# e.coli gene classification based on transcriptional factor
TF_ecoli = pd.read_excel("data/other_species/ecoli_TFSet.xlsx")
ecoli_ID_mapping = pd.read_excel("data/other_species/ecoli.xlsx")
TF_ecoli_gene_list = TF_ecoli['gene'].tolist()
TF_ecoli_gene_list0 = []
for x in TF_ecoli_gene_list:
    if ', ' in x:
        s = x.split(', ')
        TF_ecoli_gene_list0 = TF_ecoli_gene_list0 + s
    else:
        TF_ecoli_gene_list0.append(x)

TF_ecoli_gene_list0 = list(set(TF_ecoli_gene_list0))
TF_ecoli_gene_list0 = [x.lower() for x in TF_ecoli_gene_list0]
ecoli_ID_mapping['gene'] = ecoli_ID_mapping['gene'].str.lower()
ecoli_ID_mapping['TF'] = None
ecoli_ID_mapping['TF'][ecoli_ID_mapping['gene'].isin(TF_ecoli_gene_list0)] = "Yes"
ecoli_ID_mapping['TF'][~ecoli_ID_mapping['gene'].isin(TF_ecoli_gene_list0)] = "No"

# classify the gene
pro_size["TF"] = singleMapping(ecoli_ID_mapping['TF'], ecoli_ID_mapping['Entry'], pro_size["Protein"])

from scipy.stats import ttest_ind
pro_g1 = pro_size[pro_size["TF"]=="Yes"]
pro_g2 = pro_size[pro_size["TF"]=="No"]
# here just remove too long or too short amino acids?
# if not using the filter, the tendency is the same
pro_g2 = pro_g2[pro_g2["pro_length"] >= min(pro_g1["pro_length"])]
pro_g2 = pro_g2[pro_g2["pro_length"] <= max(pro_g1["pro_length"])]
# combine two pandas
pro_c = pd.concat([pro_g1, pro_g2], axis=0)


#plot
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=pro_c)
plt.xlabel("TF", fontsize=15)
plt.ylabel("volume_per_kda", fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/TF_ecoli.pdf', bbox_inches='tight')





sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])
pro_c.to_excel('data/other_species/ecoli_structure_info.xlsx')
# compare the TF in the molecular weithght
sns.catplot(x="TF", y="MW", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['MW'], pro_g2['MW'])



# re-do the above analysis using the calibrated datasets
ecoli_calibrated = pd.read_excel('data/other_species/ecoli_structure_info2.xlsx')
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=ecoli_calibrated)
sns.catplot(x="TF", y="volume_per_kda2", order=["No", "Yes"], kind="box", data=ecoli_calibrated)
plt.xlabel("TF proteins?",fontsize=15)
plt.ylabel("Volume_per_kda",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('result/Volume_per_kda_of_TF_ecoli_after_calibration.pdf', bbox_inches='tight')


#sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=ecoli_calibrated)

pro_g1 = ecoli_calibrated[ecoli_calibrated["TF"]=="Yes"]
pro_g2 = ecoli_calibrated[ecoli_calibrated["TF"]=="No"]
ttest_ind(pro_g1['volume_per_kda2'], pro_g2['volume_per_kda2'])
SUM1 = pro_g1.describe()
SUM2 = pro_g2.describe()





# how is the result if we only choose structure with high quality
ecoli_calibrated2 = ecoli_calibrated[ecoli_calibrated['score'] >=75]
sns.catplot(x="TF", y="volume_per_kda", order=["No", "Yes"], kind="box", data=ecoli_calibrated2)
plt.xlabel("TF",fontsize=15)
plt.ylabel("Volume_per_kda",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=ecoli_calibrated2)
plt.xlabel("TF",fontsize=15)
plt.ylabel("Score",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

pro_g1 = ecoli_calibrated2[ecoli_calibrated2["TF"]=="Yes"]
pro_g2 = ecoli_calibrated2[ecoli_calibrated2["TF"]=="No"]
ttest_ind(pro_g1['volume_per_kda'], pro_g2['volume_per_kda'])
ttest_ind(pro_g1['score'], pro_g2['score'])









#note：when did the enrichment analysis, it seems that the TF structure in ecoli is different from yeast.
sns.displot(pro_size, x="volume_per_kda", stat="density", common_norm=False)
plt.xlim(0.75,1.25)
plt.xlabel("Volume_per_kda", fontsize=15)
plt.ylabel("Density", fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/Volume_per_kda_ecoli.pdf', bbox_inches='tight')






pro_size = pro_size.sort_values(by=['volume_per_kda'], ascending=True)
pro_size01 = pro_size.iloc[0:200,:]
gene01= ",".join(pro_size01["Protein"].to_list())
print(gene01)

pro_size02 = pro_size.iloc[4163:4363,:]
gene02= ",".join(pro_size02["Protein"].to_list())
print(gene02)



# one interesting idea is used sce formula to calculate the protein volume in other species
pro_size["calculated_volume"] = 1.06019171e-03*pro_size["MW"] - 1.10587455
# compare the predicted and calculated for e.coli
x0 = "MW"
y0 = "Total_Volume"
y1 = "calculated_volume"
pro_size = pro_size.dropna()
x = pro_size[x0].tolist()
y = pro_size[y0].tolist()
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
fig, ax = plt.subplots(1,1,figsize=(3, 3.6))
ax.scatter(x, y, c=z, s=20)
sns.lineplot(x=x0, y=y1, data=pro_size, color='orange', linewidth=2.5,linestyle='--')
plt.xlabel(x0, fontsize=15)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
plt.savefig('result/fitted_structure_volume_ecoli.pdf', bbox_inches='tight')



