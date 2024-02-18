# Note: in this module, we will analyze the predicted protein abundance from etfl
# import self function
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import numpy

# part1
p1 = pd.read_excel("data/proteomics/enzyme_per_protein_across_compartment.xlsx")
p2 = pd.read_excel("data/proteomics/used_enzyme_per_protein_across_compartment.xlsx")

all_condition = list(p2.columns)
# analyze used and unused proteins
all_condition0 = all_condition[2:]

used_enzyme_ratio = []

initial_df = p2[["compartment"]]
for xx in all_condition0:
    print(xx)
    x0 = p1[["compartment", xx]]
    y0 = p2[["compartment", xx]]
    x0.columns = ["compartment", "total_ratio"]
    y0.columns = ["compartment", "used_ratio"]
    combine_df = pd.merge(x0, y0, left_on="compartment", right_on="compartment", how="left")
    combine_df[xx] = 1 - combine_df["used_ratio"]/combine_df["total_ratio"]
    combine_df_simple = combine_df[["compartment",xx]]
    initial_df = pd.merge(initial_df, combine_df_simple, left_on="compartment", right_on="compartment", how="left")

initial_df.to_excel("result/unused_enzyme_per_protein_in_each_compartment.xlsx")


