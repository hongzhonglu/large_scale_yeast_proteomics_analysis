# import self function
from src.protein_process import *

# input the protein abundance data
protein_copy_all1 = pd.read_excel("data/proteomics/all_protein_copy.xlsx")

# # calculate the absolute protein volume or membrane area
s1, s2 = Pro3DCal(protein_copy_all1)
s1.to_excel("data/proteomics/volume_size_across_compartment.xlsx")
s2.to_excel("data/proteomics/membrane_size_across_compartment.xlsx")

# calculate the volume ratio of proteins
s2 =Pro_3D_Volume_Ratio_Cal(protein_copy=protein_copy_all1, compartment_type="organelle") # from part 3.9
s2.to_excel("data/proteomics/volume_size_ratio_across_compartment_test.xlsx")
