import matplotlib.pyplot as plt
import os


# import self function
from src.protein_process import *

# compartment info
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
compartment_corrected = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation

df = pd.DataFrame(compartment_corrected.items(), columns=['compartment', 'gene'])
df = df.explode('gene')
# save this file
df.to_excel("data/compartment_annotation_refine.xlsx")



# compare the difference
key0 = []
len0 = []
for key in compartment.keys():
    key0.append(key)
    len0.append(len(compartment[key]))
df1 = pd.DataFrame({"compartment":key0, "gene_number": len0})


key0 = []
len0 = []
for key in compartment_corrected.keys():
    key0.append(key)
    len0.append(len(compartment_corrected[key]))
df2 = pd.DataFrame({"compartment":key0, "gene_number": len0})

# combine the dataframe
compartment_compare = pd.merge(left=df1, right=df2, left_on=['compartment'], right_on=['compartment'], how='outer')
compartment_compare.columns = ["compartment", "annotation_combine", "annotation_curation"]
compartment_compare.to_excel("data/compare_compartment_annotation_with_and_without_manual_curation.xlsx")

# save the corrected compartment annotation
mapping =[]
for key, value in compartment_corrected.items():
    print(key, value)
    new0 = [key+"@"+ x for x in value]
    mapping = mapping + new0
df = pd.DataFrame({"pair": mapping})
df1 = df['pair'].str.split('@', n=1, expand=True)
df1.columns = ['compartment','gene']
df1.to_excel("data/compartment_sce_curation.xlsx")





