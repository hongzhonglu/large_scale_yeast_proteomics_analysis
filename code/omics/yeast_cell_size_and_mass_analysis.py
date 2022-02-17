# check the following correlation between:
# 1. yeast cell size, growth rate
# 2. yeast cell mass, growth rate

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import pandas as pd

# input the physiological datasets from Rosemerry
cell_size_data = pd.read_excel("data/proteomics/yeast cell size.xls")
x0='Growth rate(h-1)'
y0='Average_cell_volume(fL)'
title0 = 'result/figure/Relation_between_yeast_cell_size_and_growth_rate_1'+ '.pdf'
print(title0)
# plt.figure()
sns.lmplot(x=x0, y=y0, data=cell_size_data, lowess=True, height=4, aspect=1, hue="data_source")
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.savefig(title0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0, 45)
plt.savefig(title0, bbox_inches='tight')

title0 = 'result/figure/Relation_between_yeast_cell_size_and_growth_rate_2' + '.pdf'
print(title0)
# plt.figure()
sns.lmplot(x=x0, y=y0, data=cell_size_data, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.savefig(title0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0, 45)
plt.savefig(title0, bbox_inches='tight')


# to calculate the cell mass under different growth rate
# assume the dry content of a yeast cell is 0.35
# assume the density of a yeast cell is 1.1126 g/ml
cell_size_data["cell_mass(pg)"] = cell_size_data['Average_cell_volume(fL)'] * 0.35 * 1.1126
x0='Growth rate(h-1)'
y0='cell_mass(pg)'
title0 = 'result/figure/Relation_between_yeast_cell_mass_and_growth_rate_1'+ '.pdf'
print(title0)
# plt.figure()
sns.lmplot(x=x0, y=y0, data=cell_size_data, lowess=True, height=4, aspect=1, hue="data_source")
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.savefig(title0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.axhline(y=13, color='k', linestyle='--')
plt.ylim(0, 20)
plt.savefig(title0, bbox_inches='tight')

title0 = 'result/figure/Relation_between_yeast_cell_mass_and_growth_rate_2' + '.pdf'
print(title0)
# plt.figure()
sns.lmplot(x=x0, y=y0, data=cell_size_data, lowess=True, height=4, aspect=1)
plt.xlabel(x0, fontsize=12)
plt.ylabel(y0, fontsize=15)
plt.savefig(title0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.axhline(y=13, color='k', linestyle='--')
plt.ylim(0, 20)
plt.savefig(title0, bbox_inches='tight')

