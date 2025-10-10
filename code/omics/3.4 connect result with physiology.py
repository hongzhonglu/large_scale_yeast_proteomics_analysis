import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# input the new absolute proteomics
mass_ratio_sce_all = pd.read_excel("data/proteomics/ProMassRatio_across_compartment_combine.xlsx")

mass_ratio_sce_all = mass_ratio_sce_all.iloc[:, 1:]


mass_ratio_sce_all0 = mass_ratio_sce_all.T

mass_ratio_sce_all0.columns = mass_ratio_sce_all0.iloc[0]
mass_ratio_sce_all1 = mass_ratio_sce_all0.drop(mass_ratio_sce_all0.index[0])


mass_ratio_sce_all1.insert(0,'sampleID', list(mass_ratio_sce_all1.index))
mass_ratio_sce_all1.insert(1,'source', [None]*len(list(mass_ratio_sce_all1.index)))
mass_ratio_sce_all1.insert(2,'growth_mode', [None]*len(list(mass_ratio_sce_all1.index)))
mass_ratio_sce_all1.insert(3,'condition_detail', [None]*len(list(mass_ratio_sce_all1.index)))
mass_ratio_sce_all1.insert(4,'dilution rate (/h)', [None]*len(list(mass_ratio_sce_all1.index)))


# input the physiological dataset
physiology = pd.read_excel("data/proteomics/physiology_collection.xlsx")

mass_ratio_sce_all1["source"] = singleMapping(physiology['source'],physiology['sampleID'],mass_ratio_sce_all1['sampleID'])
mass_ratio_sce_all1["growth_mode"] = singleMapping(physiology['growth_mode'],physiology['sampleID'],mass_ratio_sce_all1['sampleID'])
mass_ratio_sce_all1["condition_detail"] = singleMapping(physiology['condition_detail'],physiology['sampleID'],mass_ratio_sce_all1['sampleID'])
mass_ratio_sce_all1["dilution rate (/h)"] = singleMapping(physiology['dilution rate (/h)'],physiology['sampleID'],mass_ratio_sce_all1['sampleID'])
mass_ratio_sce_all1.to_excel("data/proteomics/mass_ratio_sce_all1_connected_to_physiology.xlsx")



