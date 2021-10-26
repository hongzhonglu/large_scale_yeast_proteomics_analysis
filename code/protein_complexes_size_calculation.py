# this script is to calculate the protein complexes' parameters based on single protein
# 2021-10-14


import os    ##for directory
import numpy as np
import pandas as pd
import math

# Input the datasets
single_size = pd.read_excel("result/sce_protein_size_3D_structure_combine.xlsx")

# Input the protein complexes annotation from embl database
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

complex_parse["gene"] = complex_parse["subunit"].str.replace("-MONOMER", "")






