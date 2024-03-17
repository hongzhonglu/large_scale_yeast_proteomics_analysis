import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the new absolute proteomics
mass_fraction_NCB = pd.read_excel("data/proteomics/mass_fraction_NCB.xlsx")
mass_fraction_NCB = mass_fraction_NCB.iloc[:,1:]

mass_fraction_ibrahim = pd.read_excel("data/proteomics/mass_fraction_ibrahim.xlsx")
mass_fraction_ibrahim = mass_fraction_ibrahim.iloc[:,1:]

mass_fraction_cell_system_2018 = pd.read_excel("data/proteomics/mass_fraction_cell_system_2018.xlsx")
mass_fraction_cell_system_2018 = mass_fraction_cell_system_2018.iloc[:,1:]

mass_fraction_others = pd.read_excel("data/proteomics/mass_fraction_others.xlsx")
mass_fraction_others = mass_fraction_others.iloc[:,1:]

mass_fraction_all = pd.merge(left=mass_fraction_others, right=mass_fraction_ibrahim, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_all = pd.merge(left=mass_fraction_all, right=mass_fraction_NCB, left_on=['gene'], right_on=['gene'], how="outer")
mass_fraction_all = pd.merge(left=mass_fraction_all, right=mass_fraction_cell_system_2018, left_on=['gene'], right_on=['gene'], how="outer")


# Save
mass_fraction_all.to_excel("data/proteomics/mass_fraction_combine.xlsx", index=False) # the unit the mmol/gDW

