# Some in general analysis based on total protein copies, dilution rate for datasets from different samples.
# Note here we assume that the cell average weigth is 13 pg!

import matplotlib.pyplot as plt
import os
import seaborn as sns

# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# read the protein abundance files
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx") # all the current datasets

# statistical analysis of all proteomics datasets
result_df = AllProteomicsAnalysis(pro_df=protein_copy_all1)

## start the analysis
result_df0 = result_df.transpose()
result_df0["total_copy"] = result_df0["count"]*result_df0["mean"]

# add the physiology datasets
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")




##
##
##
## the following code is wrong as physiology data contain the duplications
result_df0["sampleID"] = physiology_data['sampleID'].tolist()[0:76]
result_df0["dilution rate (/h)"] = physiology_data['dilution rate (/h)'].tolist()[0:76]
result_df0["total protein content (g/gDW)"] = physiology_data['total protein content (g/gDW)'].tolist()[0:76]
result_df0["Nitrogen source"] = physiology_data['Nitrogen source'].tolist()[0:76]
result_df0["limiting nutrient"] = physiology_data['limiting nutrient'].tolist()[0:76]

# check Jianye datasets
result_Jianye = result_df0[result_df0.index.str.contains("_M")]
result_rosemary = result_df0[result_df0.index.str.contains("prot.")]
result_rosemary_NH4_limitation = result_rosemary[result_rosemary["Nitrogen source"]=="NH4"]
result_rosemary_NH4_limitation = result_rosemary_NH4_limitation [result_rosemary_NH4_limitation ["limiting nutrient"]=="N"]


result_not_rosemary = result_df0[~result_df0.index.str.contains("prot.")]





# plot the bar plot
# Set Seaborn style
sns.set_style('darkgrid')
plt.figure(figsize = (8,4))
sns.barplot(x = "sampleID", y = "total_copy", data = result_rosemary)
plt.xticks(rotation=90)
plt.xlabel("sampleID", fontsize=12)
plt.ylabel("Total protein copy per cell", fontsize=15)
plt.show()
plt.savefig("result/figure/total_protein copy rosemary.pdf", bbox_inches='tight')


sns.set_style('darkgrid')
plt.figure(figsize = (8,4))
sns.barplot(x = "sampleID", y = "total_copy", data = result_not_rosemary)
plt.xticks(rotation=90)
plt.xlabel("sampleID", fontsize=12)
plt.ylabel("Total protein copy per cell", fontsize=15)
plt.show()
plt.savefig("result/figure/total_protein copy other datasets.pdf", bbox_inches='tight')







# check the relation between the total protein content and total protein copies
x0='total protein content (g/gDW)'
y0='total_copy'
sns.lmplot(x=x0, y=y0, data=result_df0, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# annotate the dataset points
xs=result_df0[x0].tolist()
ys=result_df0[y0].tolist()
tlab=list(result_df0.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()
plt.savefig("result/figure/total protein content-total_copy all datasets.pdf")


sns.lmplot(x=x0, y=y0, data=result_Jianye, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("result/figure/total protein content-total_copy from jianye.pdf")


sns.lmplot(x=x0, y=y0, data=result_rosemary, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("result/figure/total protein content-total_copy from rosemary.pdf")





# check the relation between the total protein content and dilution rate
x0='dilution rate (/h)'
y0='total_copy'
sns.lmplot(x=x0, y=y0, data=result_df0, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# annotate the dataset points
xs=result_df0[x0].tolist()
ys=result_df0[y0].tolist()
tlab=list(result_df0.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()
plt.savefig("result/figure/dilution rate-total_copy all datasets.pdf")


sns.lmplot(x=x0, y=y0, data=result_Jianye, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("result/figure/dilution rate-total_copy from jianye.pdf")


sns.lmplot(x=x0, y=y0, data=result_rosemary, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("result/figure/dilution rate-total_copy from rosemary.pdf")


# check the relation between the total protein content and dilution rate
# only use rosemary datasets under NH4 limitation
sns.lmplot(x=x0, y=y0, data=result_rosemary_NH4_limitation, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("result/figure/dilution rate-total_copy from rosemary NH4 limitation.pdf")




