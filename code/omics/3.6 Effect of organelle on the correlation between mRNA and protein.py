# note: not find interesting things from the following analysis

import matplotlib.pyplot as plt
import os


# import self function
from src.protein_process import *

# compartment info
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
compartment_corrected = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation

df = pd.DataFrame(compartment_corrected.items(), columns=['compartment', 'gene'])
df = df.explode('gene')

# input the physiological dataset
mRNA_protein = pd.read_excel("data/proteomics/correlation_between_protein_mRNA.xlsx")

# id mapping
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
mRNA_protein['gene'] = multiMapping(id_mapping['GeneName'], id_mapping['Entry'], mRNA_protein['Accession'])
mRNA_protein.to_excel("data/proteomics/correlation_between_protein_mRNA2.xlsx")
mRNA_protein = pd.read_excel("data/proteomics/correlation_between_protein_mRNA2.xlsx")

df['GR.Pearson.r'] = singleMapping(mRNA_protein['GR.Pearson.r'], mRNA_protein['gene'], df['gene'])
df['NM.Pearson.r'] = singleMapping(mRNA_protein['NM.Pearson.r'], mRNA_protein['gene'], df['gene'])
df['all.Pearson.r'] = singleMapping(mRNA_protein['all.Pearson.r'], mRNA_protein['gene'], df['gene'])

# remove nan values
df = df[df['GR.Pearson.r'].notna()]

column_select1 = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet',
                  'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus','nucleolus']

df0 = df[df['compartment'].isin(column_select1)]
#df0 = df[df['compartment'].str.contains("mitochondri")]
df0.to_excel("data/proteomics/correlation_between_protein_mRNA_classified based on organelle.xlsx")


