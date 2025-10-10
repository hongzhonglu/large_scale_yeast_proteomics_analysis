import matplotlib.pyplot as plt
import os


# import self function
from src.model_process import *
from src.mainFunction import *
from src.protein_process import *

# a functiion to combine the different proteomics
def combineAbosluteAbundance(omics_combine_base, omics_new, remove_column="gene"):
    df_combine_auto = pd.merge(left=omics_combine_base, right=omics_new, left_on=['all_gene'], right_on=['gene'], how="left")
    # remove the duplicated
    df_combine_auto.pop(remove_column)
    return df_combine_auto


# Part 1 Collect all the data in the unit of mmol/gDW
omics_tao1 = pd.read_excel("data/proteomics/Omics_from_tao_scale.xlsx")


# input the Carl's data under max growth
omics_carl = pd.read_excel("data/proteomics/data_PNAS_2021_scale.xlsx")
omics_carl =omics_carl[['gene','mmol/gDW']]
omics_carl.columns = ['gene', 'mmol/gDW_carl']


# input the Francesca's data under max growth
omics_francesca = pd.read_excel("data/proteomics/omics_Francesca_scale.xlsx")


# input the johan's data under four conditons
# Note: this abundance is from Nature catalysis, 2022.
# omics_johan = pd.read_excel("data/proteomics/omics_johan.xlsx")
omics_johan = pd.read_excel("data/proteomics/omics_johan_all_covered.xlsx")


# input the Tyler's data
omics_Tyler = pd.read_excel("data/proteomics/Proteome_ref.xlsx", sheet_name="Sce_0.1_Tyler")
omics_Tyler = omics_Tyler[['GeneName_S288C', 'Standard(mmol/gDW)']]
omics_Tyler.columns = ['gene', 'mmol/gDW_Tyler_D0.1']


# TO-DO add new datasets
# input data from other lab
# this data is from Proteome overabundance enables respiration but limitation onsets carbon overflow
# the unit is mmol/gDW
omics_Rahul = pd.read_excel("data/proteomics/proteomics_Rahul_2020_scale.xlsx")


# combine data from different source?
# get all genes
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
# all_gene = set(omics_carl['gene'].tolist()) | set(omics_francesca['genes'].tolist()) | set(omics_tao1['gene'].tolist()) | set(omics_johan['gene'].tolist()) | set(omics_Tyler['gene'].tolist()) | set(omics_Rahul['gene'].tolist())
# all_gene = list(set(all_gene))
all_gene = sce_gene["GeneName"].tolist()
new_df = pd.DataFrame({"all_gene": all_gene})
new_df = new_df.dropna()


#df_combine = pd.merge(left=new_df, right=omics_francesca, left_on=['all_gene'], right_on=['gene'], how="left")
#df_combine = df_combine[['all_gene','Glucose_phase(mmol/gDW)', 'Diauxic_shift(mmol/gDW)', 'Ethanol_phase(mmol/gDW)']]

df_combine = combineAbosluteAbundance(new_df, omics_francesca, remove_column="gene")


#df_combine1 = pd.merge(left=df_combine, right=omics_carl, left_on=['all_gene'], right_on=['gene'], how="left")
#df_combine1 = df_combine1[['all_gene','Glucose_phase(mmol/gDW)', 'Diauxic_shift(mmol/gDW)', 'Ethanol_phase(mmol/gDW)','mmol/gDW_carl']]

df_combine1 = combineAbosluteAbundance(df_combine, omics_carl, remove_column="gene")


df_combine2 = combineAbosluteAbundance(df_combine1, omics_tao1, remove_column="gene")
df_combine3 = combineAbosluteAbundance(df_combine2, omics_Tyler, remove_column="gene")
df_combine4 = combineAbosluteAbundance(df_combine3, omics_johan, remove_column="gene")
df_combine5 = combineAbosluteAbundance(df_combine4, omics_Rahul, remove_column="gene")

df_combine5 = df_combine5.drop('Unnamed: 0', axis=1)

"""
# get the new column
all0 = df_combine5.columns
new_columns00=[]
for x in all0:
    print(x)
    if "D=" in x:
        new_columns00.append(x)

    elif "prot." in x:
        new_columns00.append(x)

    elif "mmol" in x:
        new_columns00.append(x)
    elif "all_gene" in x:
        new_columns00.append(x)
    else:
        pass

omics_combine = df_combine5[new_columns00]"""
df_combine5.to_excel("data/proteomics/omics_measured_combine.xlsx", index=False)
df_combine5.to_csv("data/proteomics/omics_measured_combine.csv", sep='\t')

# input the new absolute proteomics
omics_cell_system_2017 = pd.read_excel("data/proteomics/omics_from_cell_systems_2017_scale.xlsx")
omics_carbon_source = pd.read_excel("data/proteomics/omics_from_carbon_source_scale.xlsx")
#omics_carbon_source.pop('ref_glc_mm_rich_aerobic(mmol/gDW)_x')
omics_jianye = pd.read_excel("data/proteomics/abundance_jianye.xlsx")
omics_kate = pd.read_excel("data/proteomics/abundance_kate.xlsx")
omics_nc_tao2 = pd.read_excel("data/proteomics/Omics_from_tao_nc_scale.xlsx")


omics_combine_auto = combineAbosluteAbundance(omics_combine_base=df_combine5, omics_new=omics_cell_system_2017, remove_column="gene")
omics_combine_auto = combineAbosluteAbundance(omics_combine_base=omics_combine_auto, omics_new=omics_carbon_source, remove_column="gene")
test = omics_cell_system_2017[omics_cell_system_2017["gene"] =="YPR080W"]
test = df_combine5[df_combine5["all_gene"] =="YPR080W"]


omics_combine_auto = combineAbosluteAbundance(omics_combine_base=omics_combine_auto, omics_new=omics_jianye, remove_column="gene")
omics_combine_auto = combineAbosluteAbundance(omics_combine_base=omics_combine_auto, omics_new=omics_kate, remove_column="gene")
omics_combine_auto = combineAbosluteAbundance(omics_combine_base=omics_combine_auto, omics_new=omics_nc_tao2, remove_column="gene")

# Save
omics_combine_auto.to_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx", index=False) # the unit the mmol/gDW
omics_combine_auto.to_csv("data/proteomics/omics_measured_combine_with_more_samples.csv", index=False) # the unit the mmol/gDW


