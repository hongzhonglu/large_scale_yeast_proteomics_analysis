# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *
import seaborn as sns



# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
# input the membrane size data
# membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment.xlsx")
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_jianye_C_limitation.xlsx")
membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])




# only take selected physiology dataset
physiology_selected = physiology_data[physiology_data["kinetic"].str.contains("_M")]
# only take selected proteomics
membrane_size_selected = membrane_size_tr0[membrane_size_tr0.index.str.contains("_M")]
membrane_size_selected["sample_ID"] = list(membrane_size_selected.index)
# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=membrane_size_selected, right=physiology_selected, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation

x = combine_data["dilution rate (/h)"].tolist()
y1 = combine_data["mitochondrial outer membrane"].tolist()
y2 = combine_data["mitochondrial inner membrane"].tolist()
fig = plt.figure(figsize=(5,5))
sns.regplot(x, y1, ci=95)
sns.regplot(x, y2, ci=95)
fig.legend(labels=['outer membrane','inner membrane'])
plt.legend(loc='upper left')
plt.xlabel('dilution rate (/h)')
plt.ylabel('Occupied membrane area (µm²)')




# check the physiological parameter correlations
plt.figure()
plt.scatter(physiology_selected["dilution rate (/h)"], physiology_selected["qGlucose (mmol/gDW h)"], marker='.', label='Glucose')
plt.scatter(physiology_selected["dilution rate (/h)"], physiology_selected["qO2 (mmol/gDW h)"], marker='.',label='O2')
plt.scatter(physiology_selected["dilution rate (/h)"], physiology_selected["qCO2 (mmol/gDW h)"], marker='.',label='CO2')
plt.scatter(physiology_selected["dilution rate (/h)"], physiology_selected["qEtOH (mmol/gDW h)"], marker='.',label='ethanol')
plt.scatter(physiology_selected["dilution rate (/h)"], physiology_selected["qAce (mmol/gDW h)"], marker='.',label='acetate')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.ylim(0, 25)
plt.show()






# only take the datasets with NH4 as nitrogen under the N limitation
# here we only explore the condition with only NH4 limitation
combine_data2 = combine_data

# check the ratio of each organelle membrane size
column_select = list(combine_data2.columns)
column_select1 = [x for x in column_select if "membrane" in x]
column_select10 = [x for x in column_select1 if "component" not in x]
column_select10 = [x for x in column_select10 if "contact site" not in x]
column_select10 = [x for x in column_select10 if "raft" not in x]
column_select10 = [x for x in column_select10 if x!="membrane"]
column_select10 = [x for x in column_select10 if "network" not in x]
column_select10 = [x for x in column_select10 if "space" not in x]

# because the size of plasma membrane and fungal-type vacuole membrane is so bigger, exclude them firstly
# column_select10 = [x for x in column_select10 if x!="plasma membrane"]
# column_select10 = [x for x in column_select10 if x!="fungal-type vacuole membrane"]
column_select10 = [x for x in column_select10 if x!="mitochondrial membrane"] # remove duplicates?
column_select10 = [x for x in column_select10 if x!="vacuolar membrane"] # remove duplicates?


membrane_only = combine_data2[column_select10]

# calculate the ratio of each organelle membrane relative to total membrane
membrane_only_ratio = membrane_only.copy()
i0 = -1
for i, x in membrane_only.iterrows():
    print(i)
    i0=i0+1
    ss = list(x)
    for j in range(len(ss)):
        print(j)
        ratio = ss[j]/sum(ss)
        membrane_only_ratio.iloc[i0,j] = ratio

membrane_only_ratio["dilution rate (/h)"] = combine_data2["dilution rate (/h)"]


# plot the figure in ratio
x0 = "dilution rate (/h)"
for y0 in column_select10:
    title0 = 'result/figure/jianye_miu_ratio_' + y0 + '.pdf'
    print(title0)
    if y0 != "mitochondrial inner membrane":
        #plt.figure(figsize=(4, 4))
        sns.lmplot(x=x0, y=y0, data=membrane_only_ratio, lowess=True, height=4, aspect=1)
        #sns.lineplot(x=x0, y=y0, data=membrane_only_ratio, marker="o")
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.xlabel(x0, fontsize=12)
        plt.ylabel(y0 + " occupied ratio", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(title0, bbox_inches='tight')
    else:
        plt.figure(figsize=(4, 4))
        sns.lineplot(x=x0, y=y0, data=membrane_only_ratio,  marker="o")
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.xlabel(x0, fontsize=12)
        plt.ylabel(y0 + " occupied ratio", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(title0, bbox_inches='tight')


# relative to the plasma's protein sectional area
membrane_to_plasma = membrane_only_ratio.copy()
for y0 in column_select10:
    membrane_to_plasma[y0] = membrane_to_plasma[y0]/membrane_only_ratio["plasma membrane"]
    membrane_to_plasma[y0] = pd.to_numeric(membrane_to_plasma[y0])


x0 = "dilution rate (/h)"
for y0 in column_select10:
    title0 = 'result/figure/jianye_miu_ratio_membrane_to_plasma_' + y0 + '.pdf'
    print(title0)
    if y0 != "mitochondrial inner membrane":
        #plt.figure(figsize=(4, 4))
        sns.lmplot(x=x0, y=y0, data=membrane_to_plasma, lowess=True, height=4, aspect=1)
        #sns.lineplot(x=x0, y=y0, data=membrane_only_ratio, marker="o")
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.xlabel(x0, fontsize=12)
        plt.ylabel(y0 + " ratio per plasma", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(title0, bbox_inches='tight')
    else:
        plt.figure(figsize=(4, 4))
        sns.lineplot(x=x0, y=y0, data=membrane_to_plasma,  marker="o")
        plt.axvline(x=0.284, color='k', linestyle='--')
        plt.xlabel(x0, fontsize=12)
        plt.ylabel(y0 + " ratio per plasma", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(title0, bbox_inches='tight')



# try to put all the result together
membrane_only_ratio["sample_ID"] = combine_data2["sample_ID"]

membrane_only_ratio2 = membrane_only_ratio.transpose()
membrane_only_ratio2.columns = combine_data2["sample_ID"]
membrane_only_ratio2 = membrane_only_ratio2.iloc[0:13,:]





# plot
x0='D=0.284_M'
y0='D=0.379_M'

membrane_only_ratio2["Relative change"] = 2*(membrane_only_ratio2[y0]-membrane_only_ratio2[x0])/(membrane_only_ratio2[x0]+membrane_only_ratio2[y0])*100
membrane_only_ratio2 = membrane_only_ratio2.sort_values(by=['Relative change'], ascending=False)
membrane_only_ratio2["compartment"] = list(membrane_only_ratio2.index)

plt.figure()
sns.set_style('darkgrid')
sns.scatterplot(x=x0, y=y0, data=membrane_only_ratio2)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.plot([0, 0.45], [0, 0.45], linewidth=2)
plt.xlim(0,0.45)
plt.ylim(0,0.45)
# annotate the dataset points
xs=membrane_only_ratio2[x0].tolist()
ys=membrane_only_ratio2[y0].tolist()
tlab=list(membrane_only_ratio2.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()
plt.savefig("result/figure/organelle membrane correlation analysis for " + y0 + " vs " + x0 + ".pdf", bbox_inches='tight')


# bar plot
x00='compartment'
y00='Relative change'
plt.figure()
sns.set_style('darkgrid')
sns.barplot(x=x00, y=y00, data=membrane_only_ratio2, capsize=.2)
plt.xlabel(x00, fontsize=12)
plt.ylabel(y00 + "_" + y0 + " vs " + x0, fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/organelle membrane relative change for " + y0 + " vs " + x0 + ".pdf", bbox_inches='tight')

