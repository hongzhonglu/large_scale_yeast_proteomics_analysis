# this script is to process compartment datasets
# 2021-11-16

# here the compartment annotation is mainly from SGD and MitoMiner

import os    # for directory
import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt

from src.mainFunction import *
from src.protein_process import *


# Input the datasets from paxDB
GO_term = pd.read_excel("data/pnas.1921890117.sd01_GO_term.xlsx")
GO_term['GO-slim mapper process term'] = GO_term['GO-slim mapper process term'].str.strip()

gene_info = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
gene_info['gene_name'] = str(gene_info['gene_name'])



GO_dict = {}
for i, x in GO_term.iterrows():
    print(i)
    ss = list(x)
    name = ss[0]
    ss = ss[1:]
    mylist = [str(x) for x in ss]
    newlist = [v for v in mylist if v != 'nan']
    # get the gene OFR name base on short gene ID

    GO_dict[name] = newlist


# build GO_term dict




# replace the gene short with long name





