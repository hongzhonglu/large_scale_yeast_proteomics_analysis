# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# 2021-11-16

import sys

# import self function
from src.mainFunction import *
from src.protein_process import *


# coeffcient change based on yeast ME model paper
# this data maybe not right!
# then how to transform the absolute abundance into one cell?
# according to this paper https://www.pnas.org/content/107/3/999, we can get the density of yeast cell 1.1029 ± 0.0026 g/mL
# assume
Vgdw = 1.7 # unit is mL/gDW Ibrahim
# average volume of a cell
Vcell = 82 #µm3
Ncell = Vgdw/Vcell*1e12 #2.07e7
# set a coefficent to transform the mmol/g.Biomass into molecular/cell
# this step should be careful in the late calculation
coefficient0 = 6.022e20/Ncell # 2.904e10




# TODO
# test the new coefficient???
coefficient1 = 6.5789e9 # this coefficient is from proteomics quality check to transfer mmol/gDW into molecular/cell

# another coefficient???
coefficient2 = 7.8298e9 # this coefficient is from cell systems paper


# To calculate the structure constraint , we need change the unit g protein per biomass into  molecular per cell!
A = 1 # mmol protein/gDW
B = A*coefficient2 # molecular protein/cell
