# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
from src.protein_process import *
import seaborn as sns

# prepare a function to analyze the protein from specific organelle
def Protein_copy_organelle(protein_copy, organelle):
    """
    :param protein_copy:
    :param organelle:
    :return:
    """
    compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    gene_nucleolus = pd.read_excel("data/gene_from_nucleolus_annotations.xlsx")
    # strictly select genes from nucleolus
    gene_nucleolus = gene_nucleolus[~gene_nucleolus['Qualifier'].str.contains('is active in')]

    # creat a dataframe to save the result
    y = organelle
    # note for nucleolus, we can also use the manually curated datasets



    # run the cycle
    if y == "plasma membrane":
        genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
    elif y == "fungal-type vacuole membrane":
        genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
        genes_select = [x for x in genes_select if
                        x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
    elif y == "endosome":
        genes_select = compartment[y]
        genes_select = [x for x in genes_select if x not in [
            "YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
    elif y == "nucleolus":
        genes_select = list(set(gene_nucleolus["gene"].tolist()))  # for the test

    else:
        genes_select = compartment[y]
    organelle_pro_copy = protein_copy[protein_copy['gene'].isin(genes_select)]
    return organelle_pro_copy


# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
# only analyze the rosemary datasets under NH4 limitation?
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
protein_copy_all_rosemary = protein_copy_all1[Sample_ID_select+["gene"]]
# curate the protein copies based on newly calculated size
df_curated = calculateCurationCoefficent()
curation_coefficent = df_curated["curation_coefficent"].to_list()
for i, sid in enumerate(Sample_ID_select):
    print(i, sid, curation_coefficent[i])
    protein_copy_all_rosemary[sid] = protein_copy_all_rosemary[sid]*curation_coefficent[i]


# find protein copy for one organelle
nucleolus_pro_copy = Protein_copy_organelle(protein_copy=protein_copy_all_rosemary, organelle='nucleolus')
nucleolus_pro_copy = nucleolus_pro_copy.dropna()


# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")


# calculate the total volume of proteins
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy_rosemary.xlsx")
# input the pro structure size data
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
pro_size = pro_size[['DBID', 'locus','Total_Volume', 'section_area_new']]

protein_copy_all1["pro_volume"] = singleMapping(pro_size['Total_Volume'],pro_size['locus'],protein_copy_all1['gene'])
sample_ID = list(protein_copy_all1.columns)
sample_ID = sample_ID[1:19]
volume_list = []
for x in sample_ID:
    ss1 = protein_copy_all1[[x,"pro_volume"]]
    ss1["value"] = ss1[x]*ss1["pro_volume"]
    ss1 = ss1[~ss1["value"].isna()]
    sum0 = sum(ss1['value'])
    # change nm^3 into um^3
    total_volume_um = sum0 / 1e9
    volume_list.append(total_volume_um)
# creat a new dataframe
total_pro_volume = pd.DataFrame({"sampleID":sample_ID,"total_pro_volume":volume_list})


nucleolus_pro_copy['pro_volume'] = singleMapping(pro_size['Total_Volume'],pro_size['locus'], nucleolus_pro_copy['gene'])
nucleolus_pro_volume = nucleolus_pro_copy.copy()
for x in sample_ID:
    ss1 = nucleolus_pro_copy[[x,"pro_volume"]]
    ss1["value"] = ss1[x]*ss1["pro_volume"] / 1e9
    nucleolus_pro_volume[x] = ss1["value"]


volume_size_tr = nucleolus_pro_volume.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[18])
volume_size_tr0 = volume_size_tr0.iloc[0:18]

# only take rosemery physiology dataset
physiology_rosemery = physiology_data[physiology_data["kinetic"].str.contains("prot.")]
# only take rosemery proteomics
volume_size_rosemery = volume_size_tr0[volume_size_tr0.index.str.contains("prot.")]
volume_size_rosemery["sample_ID"] = list(volume_size_rosemery.index)
volume_size_rosemery["total_pro_volume"] = singleMapping(total_pro_volume['total_pro_volume'], total_pro_volume['sampleID'], volume_size_rosemery["sample_ID"])
column0 = list(volume_size_rosemery.columns)[0:82] + ["total_pro_volume"] # using the automaic method,this number is 119
volume_size_rosemery_ratio = volume_size_rosemery[column0]
volume_size_rosemery_ratio1 = volume_size_rosemery_ratio.copy()

column1 = list(volume_size_rosemery.columns)[0:82]
for x in column1:
    print(x)
    volume_size_rosemery_ratio1[x] = volume_size_rosemery_ratio[x] / volume_size_rosemery_ratio["total_pro_volume"]





volume_size_rosemery_ratio1['sample_ID'] = list(volume_size_rosemery_ratio1.index)

# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery_ratio1, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")

combine_data.to_excel("data/proteomics/protein_abundance_from_organelle_Rosemary_NH4_limitation_v2.xlsx")


# further analysis
# first calculate correlation coefficient
from scipy.stats import pearsonr
correlation = []
x0 = "dilution rate (/h)"
for y0 in column1:
    corr, _ = pearsonr(combine_data[x0], combine_data[y0])
    correlation.append(corr)
new_df = pd.DataFrame({"gene":column1,"coefficient":correlation})
pd_null = new_df.sort_values(by=['coefficient'], ascending=False)
pd_null["coefficient_abs"] = abs(pd_null["coefficient"])
pd_null.to_excel("data/proteomics/protein_abundance_from_organelle_coefficient.xlsx")

# bar plot
x0='gene'
y0='coefficient'
plt.figure(figsize=(8, 6))
sns.set_style('darkgrid')
sns.barplot(x=x0, y=y0, data=pd_null, capsize=.2)
plt.xlabel(x0, fontsize=6)
plt.ylabel(y0, fontsize=15)
plt.xticks(fontsize=6)
plt.yticks(fontsize=12)
plt.xticks(rotation=90)
plt.axhline(y=0, color='k', linestyle='-')
plt.show()
plt.savefig("result/figure/organelle protein volume correlation with growth.pdf", bbox_inches='tight')


# focus on the genes with high coefficient
pd_null_filter = pd_null[pd_null['coefficient'] >= 0.9]


volume_size_rosemery_ratio2 = volume_size_rosemery_ratio1[pd_null_filter['gene'].tolist() + ['sample_ID']]
# combine the physiological datasets and proteomics datasets
combine_data = pd.merge(left=volume_size_rosemery_ratio2, right=physiology_rosemery, left_on=['sample_ID'], right_on=['kinetic'], how="left")

combine_data.to_excel("data/proteomics/protein_abundance_from_organelle_Rosemary_NH4_limitation_v2_filter.xlsx")

