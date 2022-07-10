# This module is mainly used to build a pipeline to integrate structure information with models.


# import self function
from src.mainFunction import *
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns


# for Rosemary N limitation experiments
df_curated = calculateCurationCoefficent()
# calculate the surface area
total_volume = df_curated["cell_size"].tolist()
all_radius =[(3*x/(4*math.pi))**(1/3) for x in total_volume]
cell_surface_area = [4*math.pi*x**2 for x in all_radius]
df_curated["cell_surface_area"] = cell_surface_area



# next input the organelle membrane's protein sectional area and protein volume.
# input the membrane size data
membrane_size = pd.read_excel("data/proteomics/membrane_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate
membrane_size_tr = membrane_size.transpose()
membrane_size_tr0 = membrane_size_tr.rename(columns=membrane_size_tr.iloc[1])
membrane_size_tr0 = membrane_size_tr0.iloc[2:]


# input the volume size data
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx") # curated based on cell size under different growth rate
volume_size_tr = volume_size.transpose()
volume_size_tr0 = volume_size_tr.rename(columns=volume_size_tr.iloc[1])
volume_size_tr0 = volume_size_tr0.iloc[2:]


organelle_v0 = ['mitochondrion', 'nucleus', 'cytosol',
 'endoplasmic reticulum', 'lipid droplet', 'fungal-type vacuole',
 'peroxisome', 'Golgi apparatus']

organelle_m0 = ['fungal-type vacuole membrane',
 'plasma membrane',
 'mitochondrial outer membrane',
 'endoplasmic reticulum membrane',
 'mitochondrial inner membrane',
 'Golgi membrane',
 'peroxisomal membrane',
 'nuclear membrane']

# select specific columns
membrane_size_tr1 = membrane_size_tr0[organelle_m0]
membrane_size_tr1['ID'] = list(membrane_size_tr1.index)
volume_size_tr1 = volume_size_tr0[organelle_v0]
volume_size_tr1['ID'] = list(volume_size_tr1.index)

# combine three dataframe
result1 = pd.merge(membrane_size_tr1, volume_size_tr1, how="left", on=["ID"])
result2 = pd.merge(df_curated, result1, how="left", on=["ID"])
result2.to_excel('result/combine_organelle_size_with_proteome_calculation.xlsx')




