# this script is to calculate the protein ratio in dry cell weight. This ratio should be in 0.4-0.5.
# Hongzhong Lu
# 2021-11-20

import statistics

# import self function
from src.mainFunction import *
from src.protein_process import *


def proMassRatioFromCopy(protein_copy):
    """
    The method is from kcat DL paper of Feiran and Le!
    Note: after check, it shows the transformation using data from cell system, 2017. Refer to the script from unit_unify
    :param protein_copy: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    mw["MW_Kda"] = mw["MW"] / 1000
    protein_copy["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_copy["gene"])
    protein_copy = protein_copy[protein_copy["MW"].notna()]
    protein_copy = protein_copy[protein_copy["molecular/cell"].notna()]
    # procedure 1 from kcat DL paper
    protein_copy["mmol/gDW"] = protein_copy["molecular/cell"] / (6.02 * 10e20) / 13 * 10e12
    protein_copy["mg/gDW"] = protein_copy["mmol/gDW"] * protein_copy["MW"]
    # calculate the protein ratio
    protein_ratio = sum(protein_copy["mg/gDW"]) / 1000
    return protein_ratio


def proMassRatioFromBenMethod(protein_copy):
    """
    :param protein_copy: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    mw["MW_Kda"] = mw["MW"] / 1000 # kDa = g/mmol
    protein_copy["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_copy["gene"])
    protein_copy = protein_copy[protein_copy["MW"].notna()]
    protein_copy = protein_copy[protein_copy["molecular/cell"].notna()]
    cell_volume = 32.6  #32.6 # fL/cell
    dry_content = 0.3
    cell_density = 1.1126e-12 # yeast cell density under exponential growth, [g/fL] = 1e12  g/mL

    protein_copy["mmol/cell"] = protein_copy["molecular/cell"] / 6.022e+23 * 1000
    protein_copy["mmol/gDW"] = protein_copy["mmol/cell"]/cell_volume/dry_content/cell_density

    protein_copy["mg/gDW"] = protein_copy["mmol/gDW"] * protein_copy["MW"]

    # calculate the protein ratio
    protein_ratio = sum(protein_copy["mg/gDW"]) / 1000

    # get a coefficient
    coefficent1 = 1000/6.022e+23/cell_volume/dry_content/cell_density # from molecular/cell into mmol/gDW
    coefficent2 = 1/coefficent1  # from mmol/gDW into molecular/cell
    return protein_ratio


def proMassRatioAtCell(protein_copy, yeast_cell_weight=47.65):
    """
    :param protein_copy: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    # yeast_cell_weight = 47.65 # pg, this data is based  on  one literature, which maybe not right !
    # yeast_cell_weight = 13 # pg, this data should be reasonable.
    protein_copy["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_copy["gene"])
    protein_copy = protein_copy[protein_copy["MW"].notna()]
    protein_copy = protein_copy[protein_copy["molecular/cell"].notna()]
    protein_copy["mol/cell"] = protein_copy["molecular/cell"] / 6.022e+23
    protein_copy["g/cell"] = protein_copy["mol/cell"]*protein_copy["MW"]
    protein_copy["pg/cell"] = protein_copy["mol/cell"] * protein_copy["MW"]*1e12
    # calculate the protein ratio
    protein_ratio = sum(protein_copy["pg/cell"])/yeast_cell_weight
    return protein_ratio


# should further explore how to scale protein abundance to get absolute concentrations
# input data from SGD
pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')
pro_abundance = pro_abundance[["Systematic_name","Abundance_median"]]
pro_abundance.columns = ["gene", "molecular/cell"] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
print(proMassRatioFromCopy(protein_copy=pro_abundance))
print(proMassRatioFromBenMethod(protein_copy=pro_abundance))



# input data from cell system, 2018
pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
pro_abundance = pro_abundance[["Systematic Name","Mean molecules per cell","Median molecules per cell"]]
pro_abundance.columns = ["gene", "molecular/cell","median_absolute_abundance"] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
print(proMassRatioFromCopy(protein_copy=pro_abundance))
print(proMassRatioFromBenMethod(protein_copy=pro_abundance))
print(proMassRatioAtCell(protein_copy=pro_abundance, yeast_cell_weight=13)) # from this calculation, it shown that a dry yeast cell should weight at about 13 pg. The reported yeast cell weight at about 47.65 should contain the water!


# test the protein copy under exponetional phase
protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx")

sample_name = ['Mean Copy number - Glucose','Mean Copy number - Galactose','Mean Copy number - Glycerol']

protein_copy1 = protein_copy[['gene', 'Mean Copy number - Glucose']]
protein_copy1.columns = ['gene','molecular/cell']
print(proMassRatioFromCopy(protein_copy=protein_copy1))
print(proMassRatioFromBenMethod(protein_copy=protein_copy1))


# input the data from paxdb
protein_copy = pd.read_csv("data/proteomics/abundance_table.csv")
protein_copy['molecular/cell'] = protein_copy['abundance']*80
print(proMassRatioFromCopy(protein_copy=protein_copy))




