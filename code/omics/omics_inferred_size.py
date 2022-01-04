# Note: once we have the proteomics data under different condition, then we can infer the the protein size from different sources.
# Such as we can calculate the size of complexes, the size of proteins for transporting glucose, the size of proteins from each organelle

import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

protein_abundance = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")

protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx")


































# application 1 - check the glucose transporter abundance
glucose_transporter = pd.read_excel("data/glucose_transporter_abundance_check.xlsx")
glucose_transporter = glucose_transporter[glucose_transporter['Note'].isna()]
omics_glucose_transporter = protein_abundance[protein_abundance['all_gene'].isin(glucose_transporter['gene'])]
omics_glucose_transporter.to_excel("data/omics_glucose_transporter.xlsx")


glucose_transporter_copy = protein_copy[protein_copy['gene'].isin(glucose_transporter['gene'])]
glucose_transporter_copy['section_area'] = singleMapping(glucose_transporter['section_area'],glucose_transporter['gene'],glucose_transporter['gene'])
# get the section area
glucose_transporter_copy.to_excel("data/glucose_transporter_copy.xlsx")






