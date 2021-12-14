# this script is to calculate the protein ratio in dry cell weight. This ratio should be in 0.4-0.5.
# Hongzhong Lu
# 2021-11-20

import statistics

# import self function
from src.mainFunction import *
from src.protein_process import *


def proMassRatioFromCopy(protein_abundance):
    """
    The method is from kcat DL paper!
    :param protein_abundance: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    mw["MW_Kda"] = mw["MW"] / 1000
    protein_abundance["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_abundance["gene"])
    protein_abundance = protein_abundance[protein_abundance["MW"].notna()]
    # procedure 1 from kcat DL paper
    protein_abundance["mg_gDW_m1"] = protein_abundance["absolute_abundance"] * protein_abundance["MW"] / (6.02 * 10e20) / 13 * 10e12
    # calculate the protein ratio
    protein_ratio = sum(protein_abundance["mg_gDW_m1"]) / 1000
    return protein_ratio


def proMassRatioFromCopy2(protein_abundance):
    """
    :param protein_abundance: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    mw["MW_Kda"] = mw["MW"] / 1000
    protein_abundance["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_abundance["gene"])
    protein_abundance = protein_abundance[protein_abundance["MW"].notna()]
    # procedure 1 from kcat DL paper
    protein_abundance["mg_gDW_m1"] = protein_abundance["absolute_abundance"] * protein_abundance["MW"] / (0.94*1e10)
    # calculate the protein ratio
    protein_ratio = sum(protein_abundance["mg_gDW_m1"]) / 1000
    return protein_ratio



def proMassRatioFromBenMethod(protein_abundance):
    """
    :param protein_abundance: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    mw["MW_Kda"] = mw["MW"] / 1000 # kDa = g/mmol
    protein_abundance["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_abundance["gene"])
    protein_abundance = protein_abundance[protein_abundance["MW"].notna()]
    cell_volume = 32.6 # fL/cell
    dry_content = 0.3
    cell_density = 1.1126e-12 # yeast cell density under exponential growth, [g/fL] = 1e12  g/mL

    protein_abundance["mmol/cell"] = protein_abundance["absolute_abundance"] / 6.022e+23 * 1000
    protein_abundance["mmol/g"] = protein_abundance["mmol/cell"]/cell_volume/dry_content/cell_density

    protein_abundance["mg/g"] = protein_abundance["mmol/g"] * protein_abundance["MW"]

    # calculate the protein ratio
    protein_ratio = sum(protein_abundance["mg/g"]) / 1000

    # get a coefficient
    coefficent0 = 1000/6.022e+23/cell_volume/dry_content/cell_density

    return protein_ratio



def proMassRatioAtCell(protein_abundance):
    """
    :param protein_abundance: a dataframe contain the copies of protein per cell.
    :return:
    """
    mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
    mw = mw[["locus", "proteins_molecular_weight"]]
    mw.columns = ["gene name", "MW"]
    yeast_cell = 47.65 # pg
    protein_abundance["MW"] = singleMapping(mw["MW"], mw["gene name"], protein_abundance["gene"])
    protein_abundance = protein_abundance[protein_abundance["MW"].notna()]
    protein_abundance["mol/cell"] = protein_abundance["absolute_abundance"] / 6.022e+23
    protein_abundance["g/cell"] = protein_abundance["mol/cell"]*protein_abundance["MW"]
    protein_abundance["pg/cell"] = protein_abundance["mol/cell"] * protein_abundance["MW"]*1e12
    # calculate the protein ratio
    protein_ratio = sum(protein_abundance["pg/cell"])/yeast_cell
    return protein_ratio


# should further explore how to scale protein abundance to get absolute concentrations
# input data from SGD
pro_abundance = pd.read_csv("data/sce_protein_abundance_sgd.tsv", sep='\t')
pro_abundance = pro_abundance[["Systematic_name","Abundance_median"]]
pro_abundance.columns = ["gene", "absolute_abundance"] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["absolute_abundance"].notna()]
print(proMassRatioFromCopy(protein_abundance=pro_abundance))
print(proMassRatioFromCopy2(protein_abundance=pro_abundance))



# input data from cell system, 2018
pro_abundance = pd.read_excel("data/yeast_proteomics_example_cell_system_2018.xlsx")
pro_abundance = pro_abundance[["Systematic Name","Mean molecules per cell","Median molecules per cell"]]
pro_abundance.columns = ["gene", "absolute_abundance","median_absolute_abundance"] # protein abundance per cell
pro_abundance = pro_abundance[pro_abundance["absolute_abundance"].notna()]
print(proMassRatioFromCopy(protein_abundance=pro_abundance))
print(proMassRatioFromCopy2(protein_abundance=pro_abundance))
print(proMassRatioFromBenMethod(protein_abundance=pro_abundance))
