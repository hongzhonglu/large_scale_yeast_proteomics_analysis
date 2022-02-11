# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

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

cell_volume_list = [30,31,32,33,34,35,36,37,38,39,40]
coefficent_list = [calculateCoefficient(cell_volume0=x) for x in cell_volume_list]

import matplotlib.pyplot as plt
plt.figure()
plt.scatter(cell_volume_list, coefficent_list)

