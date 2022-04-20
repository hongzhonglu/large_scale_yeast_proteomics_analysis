# this module is mainly for ecModel simulation
# which model should be used?
# first compare the model difference
# version 1
from cobra.io import load_matlab_model
from cobra import Reaction, Metabolite

import sys
import pprint

pprint.pprint(sys.path)

# import self function
from src.mainFunction import *
from src.model_process import *
dir1 = "/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecYeastGEM/ecYeastGEM_batch.mat"

ecYeast = load_matlab_model(dir1)
# a whole process to add the glucose transporter sectional area as constraints
model = ecYeast.copy()

# first input the general summary in relation between the biomass and growth rate
miu = 0.42

def calculateRatioBiomass(miu):
    protein = 0.4854 * miu + 0.3888
    lipid = -0.0456 * miu + 0.0745
    RNA = 0.1525 * miu + 0.0498
    DNA = 0.005
    Pi = 0.006
    Sulphate = 0.003
    Carbon = -0.5249 * miu + 0.4449
    return protein, lipid, RNA, DNA, Pi, Sulphate, Carbon


ss = calculateRatioBiomass(miu=0.001)



