# check the following correlation between:
# 1. organelle volume across different samples
# 2. GO-term level volume across different samples

import matplotlib.pyplot as plt
from src.protein_process import *
from sklearn.metrics import r2_score

# plot function
def linearFit(df, x_name, y_name, type):
    df = df[[x_name, y_name]]
    df = df.dropna()
    x = df[x_name]
    y = df[y_name]
    x_name = x_name.split("(")[0]
    y_name = y_name.split("(")[0]
    x_name = x_name.replace("/", "_per_")
    y_name = y_name.replace("/", "_per_")
    title0 = "result/figure/" + x_name + " vs " + y_name + " in " + type + ".pdf"
    coef = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coef)
    predict = np.poly1d(coef)
    R2 = r2_score(y, predict(x))
    print(R2)
    R2 = "{:.3f}".format(R2)
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.figure()
    plt.plot(x, y, 'yo', x, poly1d_fn(x), '--k')  # '--k'=black dashed line, 'yo' = yellow circle marker
    plt.xlabel(x_name, fontsize=18)
    plt.ylabel(y_name, fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    x_max = max(x)
    y_max = max(y)
    plt.text(x_max/3, 2*y_max/3, "R2=" + str(R2), fontsize=18)
    plt.show()
    plt.savefig(title0, bbox_inches='tight')



# input the volume size data according to organelle defition
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment_Rosemary_NH4_limitation_v2.xlsx")
colname0 = list(volume_size.columns)
colname0 = [x for x in colname0 if "prot." in x]
volume_size = volume_size[colname0]
volume_size = volume_size.iloc[:,[0,3,6,9,12,15]]


linearFit(df=volume_size, x_name="prot.1", y_name="prot.19", type="organelle")
linearFit(df=volume_size, x_name="prot.7", y_name="prot.19", type="organelle")



# input the volume size data according to GO-term
volume_size = pd.read_excel("data/proteomics/volume_size_across_go_term_Rosemary.xlsx")
colname0 = list(volume_size.columns)
colname0 = [x for x in colname0 if "prot." in x]
volume_size = volume_size[colname0]
volume_size = volume_size.iloc[:,[0,3,6,9,12,15]]

# plot
linearFit(df=volume_size, x_name="prot.1", y_name="prot.19", type="GO-term")
linearFit(df=volume_size, x_name="prot.7", y_name="prot.19", type="GO-term")







###########################
# all datasets check
###########################
# input the volume size data according to organelle defition
volume_size = pd.read_excel("data/proteomics/volume_size_across_compartment.xlsx")
print(volume_size.columns)
# plot
linearFit(df=volume_size, x_name="Min_ave_aerobic(mmol/gDW)", y_name='Min_ave_anaerobic(mmol/gDW)',type="organelle")
linearFit(df=volume_size, x_name="Glucose_phase(mmol/gDW)", y_name='Ethanol_phase(mmol/gDW)', type="organelle")
linearFit(df=volume_size, x_name="Glucose_phase(mmol/gDW)", y_name='mmol/gDW_carl', type="organelle")





# input the volume size data according to organelle defition
volume_size2 = pd.read_excel("data/proteomics/volume_size_across_go_term.xlsx")
# plot
linearFit(df=volume_size2, x_name="Min_ave_aerobic(mmol/gDW)", y_name='Min_ave_anaerobic(mmol/gDW)', type="GO-term")
linearFit(df=volume_size2, x_name="Glucose_phase(mmol/gDW)", y_name='Ethanol_phase(mmol/gDW)', type="GO-term")
linearFit(df=volume_size2, x_name="Glucose_phase(mmol/gDW)", y_name='mmol/gDW_carl', type="GO-term")


