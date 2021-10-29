# this script is to calculate the organell size parameters based on single protein
# 2021-10-14


import os    ##for directory
import numpy as np
import pandas as pd
import math

# peroxisome size doi: 10.1093/femsyr/fow038
Dp = 0.2 #um
Vp = 4/3*math.pi*(Dp/2)**3

# nucleus size,
Vn1 = 2.9 #um^3, 10.1091/mbc.E06-10-0973, 2007
Vn1_Vcell_ratio = 0.07 #, 10.1091/mbc.E06-10-0973, 2007
Dn1 = (3*Vn1/(4*math.pi))**(1/3)
Dn2 = 2 #um doi: 10.1016/j.bbagrm.2012.02.011
Vn2 = 4/3*math.pi*(Dn2/2)**3 # doi: 10.1016/j.bbagrm.2012.02.011

# vacuole size,
Vv_Vcell_ratio = 0.085 #doi: 10.1016/j.bpj.2014.03.014
Vcell = 82

# mitochondrion size
Vm_Vcell_ratio = 0.047 # https://www.pnas.org/content/117/13/7524






