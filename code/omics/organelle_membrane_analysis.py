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
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation.xlsx")
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")

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
x = combine_data["qO2 (mmol/gDW h)"].tolist()
y1 = combine_data["mitochondrial outer membrane"].tolist()
y2 = combine_data["mitochondrial inner membrane"].tolist()
fig = plt.figure(figsize=(5,5))
sns.regplot(x, y1, ci=95)
sns.regplot(x, y2, ci=95)
fig.legend(labels=['outer membrane','inner membrane'])
plt.legend(loc='upper left')
plt.xlabel('qO2 (mmol/gDW h)')
plt.ylabel('Occupied membrane area (µm²)')


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




plt.figure()
plt.scatter(combine_data["mitochondrial outer membrane"],combine_data["mitochondrial inner membrane"])
plt.xlabel('mitochondrial outer membrane')
plt.ylabel('mitochondrial inner membrane')
plt.show()



# check the physiological parameter correlations
plt.figure()
plt.scatter(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["qGlucose (mmol/gDW h)"], marker='.', label='Glucose')
plt.scatter(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["qO2 (mmol/gDW h)"], marker='.',label='O2')
plt.scatter(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["qCO2 (mmol/gDW h)"], marker='.',label='CO2')
plt.scatter(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["qEtOH (mmol/gDW h)"], marker='.',label='ethanol')
plt.scatter(physiology_rosemery["dilution rate (/h)"], physiology_rosemery["qAce (mmol/gDW h)"], marker='.',label='acetate')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.ylim(0, 25)
plt.show()






# only take the datasets with NH4 as nitrogen under the N limitation
# here we only explore the condition with only NH4 limitation
combine_data2 = combine_data[combine_data["Nitrogen source"] =="NH4"]
combine_data2 = combine_data2[combine_data2["limiting nutrient"] =="N"]
combine_data2.to_excel("data/proteomics/rosemary_data_analysis.xlsx")

# plot the figures and do the statistical analysis
x = combine_data2["qO2 (mmol/gDW h)"].tolist()
y1 = combine_data2["mitochondrial outer membrane"].tolist()
y2 = combine_data2["mitochondrial inner membrane"].tolist()
fig = plt.figure(figsize=(5,5))
sns.regplot(x, y1, ci=95)
sns.regplot(x, y2, ci=95)
fig.legend(labels=['outer membrane','inner membrane'])
plt.legend(loc='upper left')
plt.xlabel('qO2 (mmol/gDW h)')
plt.ylabel('Occupied membrane area (µm²)')
plt.axvline(x=5.4, color='k', linestyle='--')


x = combine_data2["dilution rate (/h)"].tolist()
y1 = combine_data2["mitochondrial outer membrane"].tolist()
y2 = combine_data2["mitochondrial inner membrane"].tolist()
fig = plt.figure(figsize=(5,5))
sns.regplot(x, y1, ci=95)
sns.regplot(x, y2, ci=95)
fig.legend(labels=['outer membrane','inner membrane'])
plt.legend(loc='upper left')
plt.xlabel('dilution rate (/h)')
plt.ylabel('Occupied membrane area (µm²)')
plt.axvline(x=0.18, color='k', linestyle='--')


# check the physiological parameter correlations
plt.figure(figsize=(5,5))
plt.scatter(combine_data2["dilution rate (/h)"], combine_data2["qGlucose (mmol/gDW h)"], marker='.', label='Glucose')
plt.scatter(combine_data2["dilution rate (/h)"], combine_data2["qO2 (mmol/gDW h)"], marker='.',label='O2')
plt.scatter(combine_data2["dilution rate (/h)"], combine_data2["qCO2 (mmol/gDW h)"], marker='.',label='CO2')
plt.plot(combine_data2["dilution rate (/h)"], combine_data2["qEtOH (mmol/gDW h)"], marker='.',label='ethanol')
plt.xlabel('Growth rate (/h)')
plt.ylabel('rate (mmol/gDW.h)')
plt.legend(loc='upper left')
plt.ylim(0, 25)
plt.show()






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



# plot the figure in occupied area
x0 = "dilution rate (/h)"
for y0 in column_select10:
    title0 = 'result/figure/rose_miu_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2, lowess=True, height=4, aspect=1)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied area",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')

x0 = 'qO2 (mmol/gDW h)'
for y0 in column_select10:
    title0 = 'result/figure/rose_qo2_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2, lowess=True, height=4, aspect=1)
    plt.axvline(x=5.4, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied area",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')


x0 = 'total protein content (g/gDW)'
for y0 in column_select10:
    title0 = 'result/figure/rose_total_protein_content_' + y0 + '.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=combine_data2, lowess=True, height=4, aspect=1)
    plt.axvline(x=0.286, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied area",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')


# plot the figure in ratio
x0 = "dilution rate (/h)"
for y0 in column_select10:
    title0 = 'result/figure/rose_miu_ratio_' + y0 + '.pdf'
    print(title0)

    sns.lmplot(x=x0, y=y0, data=membrane_only_ratio, lowess=True, height=4, aspect=1)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " occupied ratio",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')















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



# calculate the ratio of each organelle membrane relative to total membrane
membrane_per_cell_surface = membrane_only.copy()
df_curated = calculateCurationCoefficent()
# calculate the surface area


total_volume = df_curated["cell_size"].tolist()
all_radius =[(3*x/(4*math.pi))**(1/3) for x in total_volume]
cell_surface_area = [4*math.pi*x**2 for x in all_radius]
df_curated["cell_surface_area"] = cell_surface_area

for y0 in column_select10:
    print(y0)
    membrane_per_cell_surface[y0] = membrane_per_cell_surface[y0]/df_curated["cell_surface_area"]

membrane_per_cell_surface["dilution rate (/h)"] = combine_data2["dilution rate (/h)"]
membrane_per_cell_surface["qO2 (mmol/gDW h)"] = combine_data2["qO2 (mmol/gDW h)"]



# plot the figure in ratio
x0 = "dilution rate (/h)"
for y0 in column_select10:
    title0 = 'result/figure/rose_miu_ratio_' + y0 + '_per_cell_surface.pdf'
    print(title0)

    sns.lmplot(x=x0, y=y0, data=membrane_per_cell_surface, lowess=True, height=4, aspect=1)
    plt.axvline(x=0.18, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " per cell surface",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')



x0 = 'qO2 (mmol/gDW h)'
for y0 in column_select10:
    title0 = 'result/figure/rose_qo2_' + y0 + '_per_cell_surface.pdf'
    print(title0)
    #plt.figure()
    sns.lmplot(x=x0, y=y0, data=membrane_per_cell_surface, lowess=True, height=4, aspect=1)
    plt.axvline(x=5.4, color='k', linestyle='--')
    plt.xlabel(x0,fontsize=12)
    plt.ylabel(y0 + " per cell surface",fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(title0, bbox_inches='tight')






