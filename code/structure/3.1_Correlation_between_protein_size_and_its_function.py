# Explore how the protein structure evoluted based on their function
import matplotlib.pyplot as plt
import seaborn as sns


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
structure_quality_all = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
pro_size00 = pro_size[pro_size["locus"].isin(structure_quality_all["gene"])]

# give the quality score
pro_size00["score"] = singleMapping(structure_quality_all["score"],structure_quality_all["gene"], pro_size00["locus"])
pro_size00["MW"] = pro_size00["MW"]/1000 # change the unit as kda
pro_size00["volume_per_kda"] = pro_size00["Total_Volume"]/pro_size00["MW"]
a, b = linearFit(df=pro_size00, x_name="MW", y_name="Total_Volume")
pro_size00 = pro_size00[pro_size00["volume_per_kda"] <=2 ]
a, b = linearFit(df=pro_size00, x_name="score", y_name="volume_per_kda")
a, b = linearFit(df=pro_size00, x_name="score", y_name="pro_length")
a, b = linearFit(df=pro_size00, x_name="score", y_name="MW")
a, b = linearFit(df=pro_size00, x_name="score", y_name="Total_Volume")






# import Transcription factor
from scipy.stats import ttest_ind
TRN_sce = pd.read_excel("data/transcriptional_network/TRN_sce.xlsx")
pro_size00["TF"] = None
pro_size00["TF"][pro_size00["locus"].isin(TRN_sce["TF"])] = "Yes"
pro_size00["TF"][~pro_size00["locus"].isin(TRN_sce["TF"])] = "No"
pro_size00["id"] = singleMapping(structure_quality_all["id"],structure_quality_all["gene"],pro_size00["locus"])
pro_size00.to_excel("data/sce_protein_with_TF_classification.xlsx")


# import the normarized datasets
pro_size00 = pd.read_excel("data/sce_protein_with_TF_classification2.xlsx")
pro_g1 = pro_size00[pro_size00["TF"]=="Yes"]
pro_g2 = pro_size00[pro_size00["TF"]=="No"]
# here just remove too long or too short amino acids?
# if not using the filter, the tendency is the same
pro_g2 = pro_g2[pro_g2["pro_length"] >= min(pro_g1["pro_length"])]
pro_g2 = pro_g2[pro_g2["pro_length"] <= max(pro_g1["pro_length"])]
# combine two pandas
pro_c = pd.concat([pro_g1, pro_g2], axis=0)
sns.catplot(x="TF", y="volume_per_kda2", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['volume_per_kda2'], pro_g2['volume_per_kda2'])
sns.catplot(x="TF", y="score", order=["No", "Yes"], kind="box", data=pro_c)
SUM1 = pro_g1.describe()
SUM2 = pro_g2.describe()
# compare the TF in the molecular weithght
pro_c['MW_original'] = pro_c['MW']*1000
sns.catplot(x="TF", y="MW_original", order=["No", "Yes"], kind="box", data=pro_c)
ttest_ind(pro_g1['MW']*1000, pro_g2['MW']*1000)



# select the genes with smallest volume per kda to do enrichment analysis
pro_size00 = pro_size00.sort_values(by=['volume_per_kda'], ascending=True)
sns.displot(pro_size00, x="volume_per_kda", stat="density", common_norm=False)
plt.xlim(0.75,1.25)
pro_size01 = pro_size00.iloc[0:200,:]
gene01= ",".join(pro_size01["locus"].to_list())
print(gene01)

# select the genes with largest volume per kda to do function enrichment analysis
pro_size02 = pro_size00.iloc[5741:5941,:]
gene02= ",".join(pro_size02["locus"].to_list())
print(gene02)




# input the compartment annotation and then evaluate how the location affect the volume per kda
# compartment
compartment = getCompartmentGeneList(filter="Yes")# based on the automatic way
all_compartment = list(compartment.keys())
# select protein from nucleus
# try to input the manual check result
gene_nucleus = pd.read_excel("data/gene_belong_nucleus_annotations.xlsx")
gene_nucleus = gene_nucleus["gene"].tolist()
gene_mitochondrion = pd.read_excel("data/gene_belong_mitochondrion_annotations.xlsx")
gene_mitochondrion = gene_mitochondrion["gene"].tolist()
# here just remove the genes which belong nucleus and mitochondiron at the same time
common_gene = list(set(gene_nucleus) & set(gene_mitochondrion))
pro_size00["organelle_gene"] = "other"
pro_size00["organelle_gene"][pro_size00["locus"].isin(gene_nucleus)] = "nucleus"
pro_size00["organelle_gene"][pro_size00["locus"].isin(gene_mitochondrion)] = "mitochondrion"
pro_size00["organelle_gene"][pro_size00["locus"].isin(common_gene)] = "m & n"
sns.catplot(x="organelle_gene", y="volume_per_kda2", order=["nucleus", "mitochondrion","m & n","other"], kind="box", data=pro_size00)
plt.xlabel("Main location",fontsize=15)
plt.ylabel("Volume_per_kda",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

pro_g1 = pro_size00[pro_size00["organelle_gene"]=="mitochondrion"]
pro_g2 = pro_size00[pro_size00["organelle_gene"]=="nucleus"]
ttest_ind(pro_g1['volume_per_kda2'], pro_g2['volume_per_kda2'])




# then compare the location from nucleus
gene_nucleus = pd.read_excel("data/gene_belong_nucleus_annotations.xlsx")
gene_nucleus = gene_nucleus["gene"].tolist()
gene_TF = TRN_sce["TF"].tolist()
pro_size_n = pro_size00[pro_size00["locus"].isin(gene_nucleus)]
pro_size_n["nuclear_loc"] = "other"
pro_size_n["nuclear_loc"][pro_size_n["locus"].isin(TRN_sce["TF"])] = "TF"
sns.catplot(x="nuclear_loc", y="volume_per_kda2", order=["other","TF"], kind="box", data=pro_size_n)
plt.xlabel("Nuclear_protein_classification",fontsize=15)
plt.ylabel("Volume_per_kda",fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

pro_g1 = pro_size_n[pro_size_n["nuclear_loc"]=="TF"]
pro_g2 = pro_size_n[pro_size_n["nuclear_loc"]=="other"]
ttest_ind(pro_g1['volume_per_kda2'], pro_g2['volume_per_kda2'])





"""
normalize the datasets in R
library(readxl)
sce_TF <- read_excel("sce_protein_with_TF_classification.xlsx")
gcCount.loess <- loess(volume_per_kda~score,data=sce_TF,control = loess.control(surface = "direct"),degree=2)
predictions1<- predict (gcCount.loess, sce_TF$score)
plot(sce_TF$score,sce_TF$volume_per_kda,cex=0.1,xlab="score",ylab="volume_per_kda")
lines(sce_TF$score, predictions1,col = "red")

#sustract the influence of GC
resi <- sce_TF$volume_per_kda-predictions1
sce_TF$volume_per_kda2 <- resi
gcCount.loess <- loess(volume_per_kda2~score,data=sce_TF,control = loess.control(surface = "direct"),degree=2)
predictions2 <- predict(gcCount.loess, sce_TF$score)
plot(sce_TF$score,sce_TF$volume_per_kda2,cex=0.1,xlab="score",ylab="volume_per_kda")
lines(sce_TF$score,predictions2,col="red")
write.table(sce_TF,"sce_protein_with_TF_classification2.txt")
"""









