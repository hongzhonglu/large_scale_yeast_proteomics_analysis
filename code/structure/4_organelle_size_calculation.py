# To-do
# this script is to calculate the organell size parameters based on single protein
# 2021-10-14
import math
from src.protein_process import *

## Volume
Vcell = 82 # assume cell size is 82 um^3
Dcell = 2*(3*Vcell/(4*math.pi))**(1/3)



# peroxisome size doi: 10.1093/femsyr/fow038
Dp = 0.2 #um
Vp = 4/3*math.pi*(Dp/2)**3
Vp_Vcell_ratio = Vp/Vcell




# nucleus size,
Vn1 = 2.9 #um^3, 10.1091/mbc.E06-10-0973, 2007
Vn1_Vcell_ratio = 0.07 #, 10.1091/mbc.E06-10-0973, 2007
Dn1 = 2*(3*Vn1/(4*math.pi))**(1/3)
Dn2 = 2 #um doi: 10.1016/j.bbagrm.2012.02.011
Vn2 = 4/3*math.pi*(Dn2/2)**3 # doi: 10.1016/j.bbagrm.2012.02.011

# vacuole size,
Vv_Vcell_ratio = 0.085 #doi: 10.1016/j.bpj.2014.03.014
Vcell = 82
Vv = Vcell*0.085
Dv = 2*(3*Vv/(4*math.pi))**(1/3)

# mitochondrion size
Vm_Vcell_ratio = 0.047 # https://www.pnas.org/content/117/13/7524




## surface area
Scell = 4*math.pi*(Dcell/2)**2 # 91.27 um^2

# vacuole membrane area
# Svm_Vcell=0.2-0.3 #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4864093/pdf/nihms768320.pdf
# calculate the area based on reported ratio
Svm = Vcell*0.25  # 20.5 um^2

# calculate the area based on the above volume data
# interesting observation: the ratio of vacuole membrane area relative to the total cell surface is
# quite consistent with the reported values!!!
Svm = 4*math.pi*(Dv/2)**2 # 17.64 um^2
Svm_Scell_ratio = Svm/Scell # the ratio is 0.22, very close to the reported range 0.2-0.3


# prepare a function to calculate the ratio of surface per cell surface
# with peroxisome as an example
Vp_Vcell_ratio = 0.005/100
Vcell = 82
Vp = Vcell*Vp_Vcell_ratio
Dp = 2*(3*Vp/(4*math.pi))**(1/3)
Spm = 4*math.pi*(Dp/2)**2
Spm_Scell_ratio = Spm/Scell

getSurfaceRatio(volume_ratio = 8.5/100, Vcell = 82, Scell=91.27)
getSurfaceRatio(volume_ratio = 1/100, Vcell = 82, Scell=91.27)



# check the ratio of surface and volume in cell and protein level
## Volume
#Vcell = 40 # assume cell size is 82 um^3
#Dcell = 2*(3*Vcell/(4*math.pi))**(1/3)
Dcell = 4.24
Scell = 4*math.pi*(Dcell/2)**2 # 91.27 um^2
Vcell = 4/3*math.pi*(Dcell/2)**3
ratio_surface_to_volume1 = Scell/Vcell


# the average protein diameter 4.46 nm
Dcell = Dcell/1000
Scell = 4*math.pi*(Dcell/2)**2
S_section = math.pi*(Dcell/2)**2
Vcell = 4/3*math.pi*(Dcell/2)**3
ratio_surface_to_volume2 = Scell/Vcell
ratio_section_area_to_volume2 = S_section/Vcell