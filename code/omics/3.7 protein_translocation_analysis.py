# Note:
# all these analysis is based on the protein copy?
import seaborn as sns
import matplotlib.pyplot as plt
# import self function
from src.model_process import *
from src.protein_process import *

# input the protein mass fraction
protein_mass_fraction = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")
protein_mass_fraction = protein_mass_fraction[["gene","Min_aerobic_1","Min_aerobic_2"]]
protein_mass_fraction = protein_mass_fraction.dropna()


def create_mean_dataframe(df):
    """
    计算第2-4列每行均值，保留第一列创建新DataFrame

    参数:
        df: 原始DataFrame（至少包含4列）

    返回:
        新DataFrame: 第一列 + 2-4列均值
    """
    # 提取并保留第一列
    first_col = df.iloc[:, 0]

    # 计算第2-4列（索引1-3）每行均值
    row_means = df.iloc[:, 1:3].mean(axis=1)  # axis=1 表示按行计算[1,5](@ref)

    # 创建新DataFrame
    new_df = pd.DataFrame({
        df.columns[0]: first_col,  # 第一列数据
        'Mean_2-3': row_means  # 均值列
    })

    return new_df


# 执行计算
protein_mass_fraction0 = create_mean_dataframe(protein_mass_fraction)
protein_mass_fraction0.columns = ["gene","mass_fraction"]



# test one
protein_list = pd.read_excel("data/protein_translocation_under_different_cell_cycle.xlsx")
protein_list = protein_list.dropna()
protein_list['two_loc'] = protein_list['From'] + "_" + protein_list['To'] + "_" + protein_list['phase'].astype(str)

protein_list = protein_list[["two_loc","ORF"]]
result_dict = protein_list.groupby('two_loc')['ORF'].apply(list).to_dict()

fraction_list = []
two_loc_list = []
for key, value in result_dict.items():
    print(key, value)
    gene0 = value
    fraction0 = protein_mass_fraction0[protein_mass_fraction0['gene'].isin(gene0)]
    column_sum = fraction0['mass_fraction'].sum()
    two_loc_list.append(key)
    fraction_list.append(column_sum)

pd1 = pd.DataFrame({"two_organelle":two_loc_list,"transfer_protein_fraction":fraction_list})

# plot
# Sample list
data = fraction_list
# Create density plot
plt.figure()
sns.histplot(data)
plt.xlabel("Mass fraction")
plt.ylabel("Count")
plt.show()



# test two
# compartment info
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
compartment_corrected = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
# select the main organelle for the further analysis
organelle = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum','endosome','lipid droplet', 'fungal-type vacuole','peroxisome','ribosome','Golgi apparatus', 'plasma membrane']
# get the combination
import itertools
combinations = list(itertools.combinations(organelle, 2))
result_dict2 = {}
for x in combinations:
    print(x)
    x1 = x[0]
    x2 = x[1]
    gene_List1 = compartment_corrected[x1]
    gene_List2 = compartment_corrected[x2]
    gene_common = list(set(gene_List1) & set( gene_List2 ))
    result_dict2[x] = gene_common



fraction_list2 = []
two_loc_list2 = []
gene_number2 = []
for key, value in result_dict2.items():
    print(key, value)
    gene0 = value
    fraction0 = protein_mass_fraction0[protein_mass_fraction0['gene'].isin(gene0)]
    column_sum = fraction0['mass_fraction'].sum()
    share_protein_num = len(gene0)
    two_loc_list2.append(key)
    fraction_list2.append(column_sum)
    gene_number2.append(share_protein_num)


pd2 = pd.DataFrame({"two_organelle":two_loc_list2,"share_protein_fraction":fraction_list2, "share_protein_num":gene_number2})

# plot
# Sample list
data = fraction_list2
# Create density plot
plt.figure()
sns.histplot(data)
plt.xlabel("Mass fraction")
plt.ylabel("Count")
plt.show()
pd2.to_excel("data/organelle_interaction_data.xlsx")





# plot network graph
import networkx as nx
import matplotlib.pyplot as plt

plt.figure()
# Create an empty graph (use nx.DiGraph() for directed edges)
G = nx.Graph()
nodes = organelle
G.add_nodes_from(nodes)
edge0 = []
for m,n in pd2.iterrows():
    print(m,n)
    s0 = (n[0][0],n[0][1], n[1])
    edge0.append(s0)

weighted_edges = edge0
G.add_weighted_edges_from(weighted_edges)
weights = [G[u][v]['weight'] for u, v in G.edges()]
edge_widths = [w * 50 for w in weights]  # Adjust multiplier as needed
# Choose a layout algorithm (e.g., spring_layout, circular_layout)
pos = nx.spring_layout(G, seed=42)  # seed for reproducibility

# Draw nodes and labels
nx.draw_networkx_nodes(G, pos, node_size=200, node_color="skyblue")
nx.draw_networkx_labels(G, pos, font_size=10)

# Draw edges with proportional widths
nx.draw_networkx_edges(
    G, pos,
    width=edge_widths,  # Widths mapped to weights
    edge_color="gray",
    alpha=0.7)

#edge_labels = nx.get_edge_attributes(G, "weight")
#nx.draw_networkx_edge_labels(
#    G, pos,
#    edge_labels=edge_labels,
#    font_size=10
#)
plt.axis("off")  # Hide axes
plt.title("Weighted Network Visualization", fontsize=14)
plt.show()  #
