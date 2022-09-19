# This module is mainly used to build a pipeline to integrate structure information with models.


# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns


# second ecYeast based om deep learning
dir2 = "data/ecGEMs_and_predicted_kcat/emodel_Saccharomyces_cerevisiae_Posterior_mean.xml"
ecYeast = read_sbml_model(dir2)
gem_rxn_nov = produceRxnList(ecYeast)
gene_prot = gem_rxn_nov[gem_rxn_nov["name"].str.contains("prot_")]
gene_prot['geneID'] = gene_prot['rxnID'].str.replace("prot_", "")


#then get the protein volume information
#input the protein volume datasets
pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
gene_prot["Volume"] = singleMapping(pro_size['Total_Volume'], pro_size['locus'], gene_prot['geneID'])
gene_prot["section_area"] = singleMapping(pro_size['section_area_new'], pro_size['locus'], gene_prot['geneID'])



# plot some density graph
sns.displot(gene_prot, x="Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("Volume of single protein (nm^3)", fontsize=15)


sns.displot(gene_prot, x="section_area")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("Sectional area of single protein (nm^2)", fontsize=15)





# input the protein information in organelle level calculated from proteomics
volume_size = pd.read_excel("data/proteomics/ecGEM_volume_size_across_compartment.xlsx")
membrane_size = pd.read_excel("data/proteomics/ecGEM_membrane_size_across_compartment.xlsx")
# refine-remove some used organelles
organelle_v = collectOrganelleTerm(type="volume")
volume_size = volume_size[volume_size["compartment"].isin(organelle_v)]
volume_size = volume_size.sort_values(by=['mmol/gDW_carl'], ascending=False)
volume_size = volume_size.drop('Unnamed: 0', axis=1)
volume_size_t = volume_size.transpose()
volume_size_t.columns = volume_size_t.iloc[0]
volume_size_t = volume_size_t.iloc[1:,:]
volume_size_t = volume_size_t.apply(pd.to_numeric, errors='ignore')
# plot some density graph
organelle_v0 = ['mitochondrion', 'nucleus', 'cytosol',
 'endoplasmic reticulum', 'lipid droplet', 'fungal-type vacuole',
 'peroxisome', 'Golgi apparatus']
for xx in organelle_v0: # loop the organelle name
    sns.displot(volume_size_t, x=xx)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(xx + " protein volume (μm^3)", fontsize=15)


organelle_m = collectOrganelleTerm(type="m")
membrane_size = membrane_size[membrane_size["compartment"].isin(organelle_m)]
membrane_size = membrane_size.sort_values(by=['mmol/gDW_carl'], ascending=False)
membrane_size = membrane_size.drop('Unnamed: 0', axis=1)
membrane_size_t = membrane_size.transpose()
membrane_size_t.columns = membrane_size_t.iloc[0]
membrane_size_t = membrane_size_t.iloc[1:,:]
membrane_size_t = membrane_size_t.apply(pd.to_numeric, errors='ignore')


# plot some density graph
organelle_m0 = ['fungal-type vacuole membrane',
 'plasma membrane',
 'mitochondrial outer membrane',
 'endoplasmic reticulum membrane',
 'mitochondrial inner membrane',
 'Golgi membrane',
 'peroxisomal membrane',
 'nuclear membrane']

for xx in organelle_m0: # loop the organelle name
    sns.displot(membrane_size_t, x=xx)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(xx + " protein surface area (μm^2)", fontsize=15)






# further input the absolute protein abundance from each organelle
absolute_abundance_organelle = pd.read_excel("data/proteomics/ecGEM_absolute_pro_across_compartment.xlsx")
absolute_abundance_organelle = absolute_abundance_organelle[absolute_abundance_organelle["compartment"].isin(organelle_m0 + organelle_v0)]
absolute_abundance_organelle = absolute_abundance_organelle.sort_values(by=['mmol/gDW_carl'], ascending=False)
absolute_abundance_organelle_t = absolute_abundance_organelle.transpose()
absolute_abundance_organelle_t.columns = absolute_abundance_organelle_t.iloc[0]
absolute_abundance_organelle_t = absolute_abundance_organelle_t.iloc[1:,:]
absolute_abundance_organelle_t = absolute_abundance_organelle_t.apply(pd.to_numeric, errors='ignore')


# statistical analysis of protein abundance in organelle levels
organelle_pro_range = absolute_abundance_organelle_t.describe()
organelle_pro_range.to_excel("result/organelle_protein_abundance_range.xlsx")
absolute_abundance_organelle_t.to_excel("result/organelle_protein_abundance_ecGEMs.xlsx")



# firstly only use Rosemary datasets under N limitation
# integrate the growth phenotype datasets
Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
absolute_abundance_organelle_Rosemary = absolute_abundance_organelle_t[absolute_abundance_organelle_t.index.isin(Sample_ID_select)]
absolute_abundance_organelle_Rosemary['sampleID'] = list(absolute_abundance_organelle_Rosemary.index)
# input the physiology
# input the physiological datasets from Rosemerry
physiology_data = pd.read_excel("data/proteomics/physiology_collection.xlsx")
absolute_abundance_organelle_Rosemary['growth'] = singleMapping(physiology_data['dilution rate (/h)'],physiology_data['kinetic'],absolute_abundance_organelle_Rosemary['sampleID'] )
absolute_abundance_organelle_Rosemary['total_protein (g/gDW)'] = singleMapping(physiology_data['total protein content (g/gDW)'], physiology_data['kinetic'], absolute_abundance_organelle_Rosemary['sampleID'])
absolute_abundance_organelle_Rosemary = absolute_abundance_organelle_Rosemary.sort_values(by=['growth'], ascending=True)
pro_Rosemary = absolute_abundance_organelle_Rosemary.describe()
pro_Rosemary.to_excel("result/organelle_protein_abundance_range_rosemary.xlsx")
absolute_abundance_organelle_Rosemary.to_excel("result/organelle_protein_abundance_ecGEMs_rosemary.xlsx")





# density plot
for xx in organelle_m0: # loop the organelle name
    sns.displot(absolute_abundance_organelle_t, x=xx)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel(xx + " abs_pro abundance (mmol/gDW)", fontsize=15)
# box plot
for xx in organelle_m0 + organelle_v0: # loop the organelle name
    plt.figure(figsize=[4, 4])
    sns.set(style="darkgrid")
    sns.boxplot(y=absolute_abundance_organelle_t[xx])
    plt.xlabel(xx, fontsize=15)
    plt.ylabel("abs_pro abundance (mmol/gDW)", fontsize=15)
    plt.show()


















# Note: the following scripts were mainly used for the quality check!
# find gene based with compartment as input to compare
# the protein size, abundance within this compartment
# all metabolic genes from ecGEMs
gene_metabolic = gene_prot["geneID"].tolist()
compartment_in = organelle_v + organelle_m
m_gene_in_organelle = FingGenesForOrganelle(gene_set=gene_metabolic, compartment_list=compartment_in, compartment_type="organelle")
print(','.join(m_gene_in_organelle['plasma membrane']))
print(','.join(m_gene_in_organelle['nucleolus']))
for x in m_gene_in_organelle.keys():
    if len(m_gene_in_organelle[x]) <= 3:
        print(x)


# try to put the plasma membrane constraint into the model?
gene_select1 = m_gene_in_organelle['plasma membrane']
# get the structure based parameters
gene_prot_select1 = gene_prot[gene_prot["geneID"].isin(gene_select1)]
gene_prot_select1 = gene_prot_select1.sort_values(by=['section_area'], ascending=False)
# check the abundance for some outlier samples
#protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")
#protein_copy_all_select = protein_copy_all1[protein_copy_all1["gene"].isin(gene_select1)]



# plot some density graph
sns.displot(gene_prot_select1, x="Volume")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("Volume of single protein (nm^3)", fontsize=15)


sns.displot(gene_prot_select1, x="section_area")
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel("Sectional area of single protein (nm^2)", fontsize=15)


