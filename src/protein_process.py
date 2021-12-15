import math
import numpy as np
import statistics
import pandas as pd
from src.mainFunction import *

def getSurfaceRatio(volume_ratio = 0.005/100, Vcell = 82, Scell=91.27):
    """
    The function is used to calculate the surface area ratio between organelle and cell.
    :param volume_ratio:
    :param Vcell:
    :param Scell:
    :return:
    """
    Vany = Vcell * volume_ratio
    Dany = 2 * (3 * Vany / (4 * math.pi)) ** (1 / 3)
    Sanym = 4 * math.pi * (Dany / 2) ** 2
    ratio = Sanym / Scell
    return ratio

def getGeneListFromLocation(gene_location_annotation, location):
    """
    The function is to extract gene list based on its compartment information
    :param gene_location_annotation:
    :param location:
    :return:
    """
    #location = 'mitochondrial envelope'
    gene_subset = gene_location_annotation[gene_location_annotation["GO_Name"]==location]
    gene_list = list(set(gene_subset["Systematic_name"].tolist()))
    return gene_list

def getStructureSize(genes_select0, pro_size0, pro_abundance0):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])
    combine_df["abundance"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],
                                            combine_df["gene"])

    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    abundance0 = combine_df["abundance"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    abundance_median = statistics.median(abundance1) # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
    abundance_update = []
    for x in abundance0:
        if np.isnan(x) == False:
            x0 = x
        else:
            x0 = abundance_median
        abundance_update.append(x0)
    combine_df["abundance_update"] = abundance_update

    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["abundance_update"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["abundance_update"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    return total_volume_um, total_area_um


def splitAbundance(pro_df):
    """
    The function is used to quality check of protein abundance in molecular/cell before entering next step.
    :param pro_df: A dataframe should columns-gene,molecular/cell.
    :return:
    """
    # sometimes it shows mutiple proteins together have one abundance value, here we need a function to do the quality check!
    len1 = pro_df.shape[0]
    pro_df1 = pro_df[pro_df['gene'].str.contains(';')]
    if (len(pro_df1) > 0):
        print('Mutiple protein have one abudance value! Need quality check.')
    pro_df2 = pro_df[~pro_df['gene'].str.contains(';')]
    gene0 = []
    abundance0 = []
    for i, x in pro_df1.iterrows():
        # print(i, x)
        s = x["gene"].split(";")
        len0 = len(s)
        gene0 = gene0 + s
        v = [x["molecular/cell"] / len0] * len0
        abundance0 = abundance0 + v
    pro_df1 = pd.DataFrame({"gene": gene0, "molecular/cell": abundance0})

    pro_df = pd.concat([pro_df1, pro_df2], axis=0)
    len2 = pro_df.shape[0]

    if (len2 > len1):
        print('Complete the quality check!')

    return pro_df




