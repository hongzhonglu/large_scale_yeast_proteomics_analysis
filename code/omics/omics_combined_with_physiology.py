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
membrane_size = pd.read_excel("data/proteomics/membrance_size_across_compartment.xlsx")
membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])




# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
membrane_size_rosemery = membrane_size_tr0[membrane_size_tr0.index.str.contains("prot.")]
membrane_size_rosemery["sample_ID"] = list(membrane_size_rosemery.index)
# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=membrane_size_rosemery, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")
# further filter based on Nitrogen limitation or carbon limitation


# plot the figures and do the statistical analysis
plt.figure()
plt.scatter(combine_data["qO2 (mmol/gDW h)"],combine_data["mitochondrial outer membrane"])
plt.xlabel('qO2 (mmol/gDW h)')
plt.ylabel('mitochondrial outer membrane')
plt.show()

plt.figure()
plt.scatter(combine_data["qO2 (mmol/gDW h)"], combine_data["mitochondrial inner membrane"])
plt.xlabel('qO2 (mmol/gDW h)')
plt.ylabel('mitochondrial inner membrane')
plt.show()


plt.figure()
plt.scatter(combine_data["mitochondrial outer membrane"],combine_data["mitochondrial inner membrane"])
plt.xlabel('mitochondrial outer membrane')
plt.ylabel('mitochondrial inner membrane')
plt.show()




x = combine_data["qO2 (mmol/gDW h)"].tolist()
y1 = combine_data["mitochondrial outer membrane"].tolist()
y2 = combine_data["mitochondrial inner membrane"].tolist()
fig = plt.figure(figsize=(4,6))
sns.regplot(x, y1, ci=95)
sns.regplot(x, y2, ci=95)
fig.legend(labels=['outer membrane','inner membrane'])
plt.legend(loc='upper left')
plt.xlabel('qO2 (mmol/gDW h)')
plt.ylabel('Occupied membrane area (µm²)')






# further explore the correlation between inner membrane and outer membrane
membrane_size_tr0 = membrane_size_tr0.iloc[2:]
# if only take datasets from sysbio in recent years
physiology_sysbio = physiology_data[physiology_data["source"]=="sysbio"]


sample_sysbio = physiology_sysbio["kinetic"].tolist()
membrane_size_sysbio = membrane_size_tr0[membrane_size_tr0.index.isin(sample_sysbio)]
# further remove Jianye datasets
sample_sysbio2 = [x for x in sample_sysbio if "D=" not in x]
membrane_size_sysbio = membrane_size_sysbio[membrane_size_sysbio.index.isin(sample_sysbio2)]
membrane_size_sysbio0 = membrane_size_sysbio[["mitochondrial outer membrane", "mitochondrial inner membrane"]]
membrane_size_sysbio0.to_excel("data/proteomics/membrane_size_sysbio.xlsx")


# plot
plt.figure()
plt.scatter(membrane_size_sysbio["mitochondrial outer membrane"], membrane_size_sysbio["mitochondrial inner membrane"])
plt.xlabel('mitochondrial outer membrane')
plt.ylabel('mitochondrial inner membrane')


xs=membrane_size_sysbio["mitochondrial outer membrane"].tolist()
ys=membrane_size_sysbio["mitochondrial inner membrane"].tolist()
tlab=list(membrane_size_sysbio.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()


fig = plt.figure(figsize=(4,6))
sns.regplot(xs, ys, ci=95)
plt.xlabel('mitochondrial outer membrane')
plt.ylabel('mitochondrial inner membrane')


ax = sns.regplot(xs, ys, ci=95)
ax.set(yscale='log')

# fit a nonparametric regression using a lowess smoother.
sns.lmplot(x="mitochondrial outer membrane", y="mitochondrial inner membrane", data=membrane_size_sysbio,
           lowess=True)



