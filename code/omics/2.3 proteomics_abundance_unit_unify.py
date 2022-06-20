# The script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# based on the cellular volume and density
# In general, it could assume that a cell volume is 32.6 fL (or a cell weigth is 13 pg).
# 2021-11-16


import pandas as pd


# coeffcient change based on yeast ME model paper
# this data maybe not right!
# then how to transform the absolute abundance into one cell?
# according to this paper https://www.pnas.org/content/107/3/999, we can get the density of yeast cell 1.1029 ± 0.0026 g/mL
# assume
Vgdw = 1.7 # unit is mL/gDW based on Ibrahim and Metabolic Engineering 13 (2011) 294–306
# average volume of a cell
Vcell = 82 #82 #µm3
Ncell = Vgdw/Vcell*1e12 #2.07e7
# set a coefficent to transform the mmol/g.Biomass into molecular/cell
# this step should be careful in the late calculation
coefficient0 = 6.022e20/Ncell #2.904e10




# TODO
# test the new coefficient???
coefficient1 = 6.5789e9 # based on Ben. this coefficient is from proteomics quality check to transfer mmol/gDW into molecular/cell

# another coefficient???
coefficient2 = 7.8298e9 # this coefficient is from cell systems paper


# To calculate the structure constraint , we need change the unit g protein per biomass into  molecular per cell!
A = 1 # mmol protein/gDW
B = A*coefficient2 # molecular protein/cell




# in-depth analysis based on Ben
cell_volume = 32.6  # 32.6 # fL/cell
dry_content = 0.35  # https://onlinelibrary.wiley.com/doi/pdf/10.1002/j.2050-0416.1952.tb02660.x#:~:text=Yeast%20cakes%20produced%20by%20normal,the%20conditions%20of%20growth%2C%20to
cell_density = 1.1126e-12  # g/fL yeast cell density under exponential growth, [g/fL] = 1e12  g/mL
coefficent10 = 1000 / 6.022e+23 / cell_volume / dry_content / cell_density  # from molecular/cell into mmol/gDW
coefficent20 = 1 / coefficent10  # from mmol/gDW into molecular/cell


def calculateCoefficient(cell_volume0):
    cell_volume = cell_volume0  # 32.6 # fL/cell
    dry_content = 0.35  # https://onlinelibrary.wiley.com/doi/pdf/10.1002/j.2050-0416.1952.tb02660.x#:~:text=Yeast%20cakes%20produced%20by%20normal,the%20conditions%20of%20growth%2C%20to
    cell_density = 1.1126e-12  # g/fL yeast cell density under exponential growth, [g/fL] = 1e12  g/mL
    coefficent10 = 1000 / 6.022e+23 / cell_volume / dry_content / cell_density  # from molecular/cell into mmol/gDW
    coefficent20 = 1 / coefficent10  # from mmol/gDW into molecular/cell
    return coefficent20

def calculateCurationCoefficent():
    # curation of rosemary datasets based on the fitted cell volume under different growth rates
    # fitting formula to calculate the coefficients
    # when miu > 0.2, cell_volume = 77.32 miu + 15.771
    # when miu < 0.2, average volume is 28 um^3

    # Rosemary sample ID
    Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
    coefficient2 = 7.8298e9  # this is the original coefficient used in cell system paper!
    growth_rate = [0.05, 0.1, 0.13, 0.18, 0.3, 0.35]
    cell_volume_fit = [28, 28, 28, 28, 38.967, 42.833]
    coefficent_list = [calculateCoefficient(cell_volume0=x) for x in cell_volume_fit]
    curation_coefficent = [x / coefficient2 for x in coefficent_list]
    new_coefficent = []
    for x in curation_coefficent:
        print(x)
        s = [x] * 3
        new_coefficent = new_coefficent + s
    cell_volume_all = []
    for x in cell_volume_fit:
        print(x)
        s = [x] * 3
        cell_volume_all = cell_volume_all + s
    curation_info_rosemary = pd.DataFrame({"ID": Sample_ID_select, "cell_size": cell_volume_all, "curation_coefficent": new_coefficent})

    return curation_info_rosemary

# test
df_curated = calculateCurationCoefficent()
