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


def setReferenceProCopy():
    """
    This function is used to set the reference protein molecular/cell! If we still could not find the values from the
    reference value, then we can use the values at 5 percentile or 10 percentile?
    :return:
    """
    pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
    pro_abundance = pro_abundance[["Systematic Name", "Mean molecules per cell", "Median molecules per cell"]]
    pro_abundance.columns = ["gene", "absolute_abundance", "median_absolute_abundance"]  # protein abundance per cell
    pro_abundance.columns = ['gene', 'molecular/cell', "median_absolute_abundance"]
    pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
    reference_copy = pro_abundance
    statistics_analysis = reference_copy.describe()
    v_five_percent = reference_copy['molecular/cell'].quantile(0.05)
    v_ten_percent = reference_copy['molecular/cell'].quantile(0.1)
    return reference_copy, v_five_percent, v_ten_percent


def getProAundance(genes_select0, pro_abundance0):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    combine_df["molecular/cell"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],
                                            combine_df["gene"])
    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    abundance0 = combine_df["molecular/cell"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    abundance_median = statistics.median(abundance1) # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
    abundance_update = []
    for x in abundance0:
        if np.isnan(x) == False:
            x0 = x
        else:
            x0 = abundance_median
        abundance_update.append(x0)
    combine_df["molecular/cell_local"] = abundance_update
    
    # use the second choice
    # load the reference molecular copies
    ref_abundance, v_5, v_10 = setReferenceProCopy()
    # set a dict
    gene_abundance = {}
    for i, x in ref_abundance.iterrows():
        print(i)
        gene_abundance[x['gene']] = x['molecular/cell']
    
    abundance_update2 = []
    for i, x in combine_df.iterrows():
        print(i)
        ss = x['molecular/cell']
        if np.isnan(ss) == False:
            x0 = ss
        elif x['gene'] in ref_abundance['gene'].tolist():
            x0 = gene_abundance[x['gene']]   
        else:
            x0 = v_5
        abundance_update2.append(x0)
    
    combine_df["molecular/cell_global"] = abundance_update2

    return combine_df


def getStructureSize(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0
    combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])
    combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])

    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell_global"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell_global"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_volume_um, total_area_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell_global"] * combine_df["Volume"]
        combine_df["total_area"] = combine_df["molecular/cell_global"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_volume_um, total_area_um, combine_df



def splitAbundance(pro_df):
    """
    The function is used to quality check of protein abundance in molecular/cell or other unit before entering next step.
    :param pro_df: A dataframe should columns-gene,molecular/cell.
    :return:
    """
    # sometimes it shows mutiple proteins together have one abundance value, here we need a function to do the quality check!
    len1 = pro_df.shape[0]
    colnames=pro_df.columns
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
        v = [x[colnames[1]] / len0] * len0
        abundance0 = abundance0 + v
    gene0=[x.strip(" ") for x in gene0]
    pro_df1 = pd.DataFrame({"gene": gene0, colnames[1]: abundance0})

    pro_df = pd.concat([pro_df1, pro_df2], axis=0)
    len2 = pro_df.shape[0]

    if (len2 > len1):
        print('Complete the quality check!')

    return pro_df




