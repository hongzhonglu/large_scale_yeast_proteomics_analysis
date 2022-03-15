# this is data from Ben!
# 2022.03.14

import pandas as pd

# Extract of supplementary table 5:
data = pd.read_csv("data/example_data/s05_proteomics_data_raw.csv", index_col=0)

# Molecular weights:
MW = data["molecular_weight"] # Da = g/mol
MW = MW / 1000 # kDa = g/mmol
print(MW)

# Rest of the data (abundances + uncertainties):
data = data.iloc[:,1:]
print(data)


# Cell volumes (fL/cell) (extract of supplementary table 23):
cell_volumes = pd.read_csv("data/example_data/s23_cell_volume.csv", index_col=0)
print(cell_volumes)

# convert unit
for (col_name, d) in data.iteritems():
    if col_name.endswith("uncertainty"):
        mean_name = col_name.replace("uncertainty", "mean")
        data[col_name] = data[col_name] / 100 * data[mean_name]

print(data)






# Convert values to mmol/cell:
data = data / 6.022e+23 * 1000

# Iterate through the dataset and divide by the corresponding cell volume, to get mmol/fL:
for (col_name, d) in data.iteritems():
    chemo_name = col_name.replace("_uncertainty", "").replace("_mean", "")
    data[col_name] = data[col_name] / cell_volumes.loc[chemo_name]["cell_volume"]

# Finally, convert to mmol/gDW:
dry_content = 0.3
cell_density = 1.105e-12
data = data / cell_density / dry_content
print(data)

# data validation
for (c, vol) in cell_volumes.iteritems():
    print(vol * cell_density * dry_content * 1e12)  # fL/cell * g/fL * gDW/g * pgDW/gDW = pgDW/cell

for (col_name, col_data) in data.iteritems():
    prot_fraction = sum(col_data * MW)  # mmol/gDW * g/mmol = g/gDW
    print(col_name + ": " + str(prot_fraction))

data.to_csv("data/example_data/s05_proteomics_data_processed.csv")





