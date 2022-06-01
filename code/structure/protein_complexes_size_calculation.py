# this script is to calculate the protein complexes' parameters based on single protein
# 2021-10-14


import os
import numpy as np
import pandas as pd
import math
from src.mainFunction import *


# Input the datasets
single_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")

# Input the protein complexes annotation from EBML database
protein_complex = pd.read_csv("data/info_yeast/Yeast_complex_portal.csv")
complex_list=[]
subunit_list=[]
subunit_count=[]
for i, x in protein_complex.iterrows():
    print(i,x["Component coefficients"])
    complex_name = x["Product Name"]
    ss = x["Component coefficients"]
    ss = ss.replace("TUPLE //", "").replace(" ", "")
    ss1 = ss.split("//")
    odd_numbers = [y for x, y in enumerate(ss1) if x % 2 != 0]
    even_numbers = [y for x, y in enumerate(ss1) if x % 2 == 0]
    complex_name2 = [complex_name]*len(even_numbers)
    complex_list = complex_list + complex_name2
    subunit_list = subunit_list + even_numbers
    subunit_count = subunit_count + odd_numbers

complex_parse = pd.DataFrame({"complex":complex_list,"subunit":subunit_list,"count":subunit_count})
complex_parse["count"] = complex_parse["count"].astype(float)
complex_parse.to_excel("data/complex_info.xlsx")


complex_parse["gene"] = complex_parse["subunit"].str.replace("-MONOMER", "")
complex_parse["Total_Volume"] = singleMapping(single_size["Total_Volume"],single_size["locus"],complex_parse["gene"])
#complex_parse["radius_new"] = singleMapping(single_size["radius"],single_size["locus"],complex_parse["gene"])
complex_parse["section_area_new"] = singleMapping(single_size["section_area_new"],single_size["locus"],complex_parse["gene"])

# calculate the protein complex volume
all_complexes = list(set(complex_parse["complex"].tolist()))
total_volume = []
# total_section_area ??
for i in all_complexes:
    print(i)
    complex0 = complex_parse[complex_parse["complex"] == i]
    coefficient = complex0["count"].tolist()
    each_volume = complex0["Total_Volume"].tolist()
    volume0 =[x*y for x, y in zip(coefficient, each_volume)]
    volume_add = sum(volume0)
    total_volume.append(volume_add)

complex_volume = pd.DataFrame({"complex": all_complexes,"Total_Volume":total_volume})
complex_parse["Volume_complex"] = singleMapping(complex_volume["Total_Volume"],complex_volume["complex"],complex_parse["complex"])
# it shows that fatty acid synthetase with biggest size
complex_parse1 = complex_parse.sort_values(by='Volume_complex')
complex_volume = complex_volume.sort_values(by='Total_Volume')
ax = complex_volume["Total_Volume"].plot.hist(bins=12, alpha=0.5)
ax.set_title("Complex protein volume distribution")
ax.set_xlabel("Complex protein volume")
ax.set_ylabel("Density")












# build connects between complexes and yeast-GEM
yeast_gem_gene = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="GENES")
complex_parse["existence_in_GEM"] = complex_parse["gene"].isin(yeast_gem_gene["NAME"])
# It seems that complex annotation from EBML database is not fully consisted with yeast-GEM
# also it shows that the part of subunit is in yeast-GEM while other is not in yeast-GEM, that means
# the annotation from EBML may be not correct.
# On the other hand, the complex annotation from EBML can be further used to update yeast-GEM
complex_parse.to_excel("result/complex_ebml_with_gem.xlsx")


