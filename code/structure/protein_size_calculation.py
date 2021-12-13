# this script is to calculate the protein parameters based on rough estimation as well as protein 3D structure
# Hongzhong Lu
# 2021-10-14


import os    ##for directory
import numpy as np
import pandas as pd
import math
from src.mainFunction import *

# Input the datasets
data_merge = pd.read_csv("data/sce_protein_weight.tsv",sep='\t')

# calculate the size of proteins
radius_list = []
volume_list = []
section_area_list = []
for i, x in data_merge.iterrows():
    print(x["proteins_molecular_weight"])
    molecular_weight = x["proteins_molecular_weight"]
    radius = 0.066*molecular_weight**(1/3)
    radius_list.append(radius) # unit is the nm!!
    volume = 4/3*math.pi*radius**3
    volume_list.append(volume)
    section_area = math.pi*radius**2
    section_area_list.append(section_area)


data_merge["radius"] = radius_list
data_merge["volume"] = volume_list
data_merge["section_area"] = section_area_list

# second part
# calculate the protein size based on its protein 3D structures
input_file = "data/OutputDir_alphafold_pdb_11554388610859884589.txt"
with open(input_file) as file_in:
    lines = []
    for line in file_in:
        lines.append(line)
lines0 = [x for x in lines if "Reading hydrogens is turned on" not in x]
lines1 = lines0[7:]
p1 = []
p2 = []
p3 = []
p4 = []
p5 = []
p6 = []
for x in lines1:
    print(x)
    ss = x.split(" ")
    ss0 = [x for x in ss if x is not '']
    ss0 = [x.replace("\n", "") for x in ss0]
    ss0 = [x.replace(",", ".") for x in ss0]
    s1 = ss0[0]
    s2 = float(ss0[1])
    s3 = float(ss0[2])
    s4 = float(ss0[3])
    s5 = float(ss0[4])
    s6 = float(ss0[5])
    p1.append(s1)
    p2.append(s2)
    p3.append(s3)
    p4.append(s4)
    p5.append(s5)
    p6.append(s6)
# note the original unit for the volume is Å
# 1Å = 0.1 nm; 1Å^3 = 0.001 nm^3
volume_df = pd.DataFrame({"Protein":p1,"Total_Volume":p2,"Void_Volume":p3,"VDW_Volume":p4,"Packing Density":p5,"Time_Taken_ms":p6})

# change the unit from Å to nm
volume_df0 = volume_df.copy()
volume_df0["Total_Volume"] = volume_df["Total_Volume"]/1000
volume_df0["Void_Volume"] = volume_df["Void_Volume"]/1000
volume_df0["VDW_Volume"] = volume_df["VDW_Volume"]/1000


# calculate the radius and sectional area of a protein
total_volume = volume_df0["Total_Volume"].tolist()
all_radius =[(3*x/(4*math.pi))**(1/3) for x in total_volume]
all_sectional_area = [math.pi*x*x for x in all_radius]
volume_df0["radius"] = all_radius
volume_df0["section_area"] = all_sectional_area


# read id map file and merge the protein volume calculated by different methods
id_mapping = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
id_mapping["id"] = id_mapping["id"].str.replace(".pdb","")

data_merge["Protein"] = singleMapping(id_mapping["id"], id_mapping["gene"], data_merge["locus"])

data_merge1 = data_merge[~data_merge['Protein'].isnull()]
data_merge1["Total_Volume"] = singleMapping(volume_df0["Total_Volume"],volume_df0["Protein"],data_merge1["Protein"])
data_merge1["Void_Volume"] = singleMapping(volume_df0["Void_Volume"],volume_df0["Protein"],data_merge1["Protein"])
data_merge1["VDW_Volume"] = singleMapping(volume_df0["VDW_Volume"],volume_df0["Protein"],data_merge1["Protein"])
data_merge1["radius_new"] = singleMapping(volume_df0["radius"],volume_df0["Protein"],data_merge1["Protein"])
data_merge1["section_area_new"] = singleMapping(volume_df0["section_area"],volume_df0["Protein"],data_merge1["Protein"])

CHECK = data_merge1[data_merge1["locus"]=="YMR231W"] # the protein has no structure size data
data_merge1=data_merge1[~(data_merge1["locus"]=="YMR231W")]



# plot the density graph
ax = data_merge1.plot.scatter(x='volume', y='Total_Volume')
ax.set_title("")
ax.set_xlabel("Rough estimated volume(nm^3)")
ax.set_ylabel("Structure_based volume(nm^3)")
ax.set_xlim(0,350)
ax.set_ylim(0,350)

# however some proteins did not have protein 3D structures, then we still need the data from the molecualr weight
data_merge1 = data_merge1[['DBID', 'locus','Total_Volume', 'section_area_new']]
gene_no_structure = list(set(data_merge["locus"])-set(data_merge1["locus"]))
data_need_added = data_merge[data_merge["locus"].isin(gene_no_structure)]
data_need_added = data_need_added[['DBID', 'locus','volume', 'section_area']]
data_need_added.columns = ['DBID', 'locus','Total_Volume', 'section_area_new']

# combine the dataset
data_for_save = pd.concat([data_merge1, data_need_added], axis=0)
data_for_save.to_excel("result/sce_protein_size_3D_structure.xlsx")



