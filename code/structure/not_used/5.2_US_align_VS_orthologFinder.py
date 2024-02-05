# orthlogFinder process
import pandas as pd
dir0 = '/Users/xluhon/Documents/data_for_structure_align/protein_seq_for_alphafold/OrthoFinder/Results_Sep02/Orthogroups/Orthogroups.xlsx'
ortholog_relation = pd.read_excel(dir0)
ortholog_relation0 = ortholog_relation.dropna()


us_align_dir0 = '/Users/xluhon/Documents/tm_out_two_strains/'

pair = []
tm_score = []
identity_us_align = []
#process:
with open("data/organism_result_two_strains.sh", "w") as outfile:
    for i, x in ortholog_relation0.iterrows():
        #print(i, x)
        y1 = x[1].split(', ')
        y2 = x[2].split(', ')

        # check the name as the protein stucture ID
        y1 = [x.split('|')[1] for x in y1]
        y2 = [x.split('|')[1] for x in y2]

        y1 = ['AF-' + x + '-F1-model_v3.pdb' for x in y1]
        y2 = ['AF-' + x + '-F1-model_v3.pdb' for x in y2]

        # try to save the file
        for i in y1:
            for j in y2:
                s_out = i + '@@' + j
                print(s_out)
                file0 = us_align_dir0 + s_out + '.txt'
                file01 = open(file0).readlines()
                if len(file01) >= 2:
                    tm_score1 = file01[14]
                    tm_score2 = file01[15]
                    identity = float(file01[13].split("= ")[3])
                    v1 = float(tm_score1.split(' (')[0].split('= ')[1])
                    v2 = float(tm_score2.split(' (')[0].split('= ')[1])
                    mean0 = (v1 + v2) / 2
                    pair.append(s_out)
                    tm_score.append(mean0)
                    identity_us_align.append(identity)

result_df = pd.DataFrame({"pair_id":pair,"tm_score": tm_score, "identity_us_align":identity_us_align})
result_df1 = result_df['pair_id'].str.split('@@', 1, expand=True)
result_df1.columns = ['pro1', 'pro2']



def find_between_r(s, first, last ):
    out = []
    for x in s:
        try:
            start = x.rindex(first) + len(first)
            end = x.rindex(last, start)
            s = x[start:end]
            out.append(s)
        except ValueError:
            s = ""
            out.append(s)
    return out




result_df2 = result_df1.copy()

result_df2['pro1'] = find_between_r(s=result_df1['pro1'].tolist(), first='AF-', last='-F1')
result_df2['pro2'] = find_between_r(s=result_df1['pro2'].tolist(), first='AF-', last='-F1')

result_df['pro1'] = result_df2['pro1']
result_df['pro2'] = result_df2['pro2']

result_df['combine1'] = result_df['pro1'] + "@@" + result_df['pro2']
result_df['combine2'] = result_df['pro2'] + "@@" + result_df['pro1']




# read blast
canda_yeast = pd.read_csv('/Users/xluhon/Documents/data_for_structure_align/seq_blast_result/canda_yeast.txt', sep="\t", header= None)
ss0 = ["geneID", "hitID", "pident", "length", "mismatch", "gapopen",
       "qstart", "qend", "sstart", "send", "evalues", "bitscore"]
canda_yeast.columns = ss0
s1 = canda_yeast["geneID"].tolist()
s10 = [x.split("|")[1] for x in s1]
canda_yeast["geneID"] = s10

s2 = canda_yeast["hitID"].tolist()
s20 = [x.split("|")[1] for x in s2]
canda_yeast["hitID"] = s20

canda_yeast["combine"] = canda_yeast["geneID"] + "@@" + canda_yeast["hitID"]






yeast_canda = pd.read_csv('/Users/xluhon/Documents/data_for_structure_align/seq_blast_result/yeast_canda.txt', sep="\t", header= None)
ss0 = ["geneID", "hitID", "pident", "length", "mismatch", "gapopen",
       "qstart", "qend", "sstart", "send", "evalues", "bitscore"]
yeast_canda.columns = ss0


s1 = yeast_canda["geneID"].tolist()
s10 = [x.split("|")[1] for x in s1]
yeast_canda["geneID"] = s10

s2 = yeast_canda["hitID"].tolist()
s20 = [x.split("|")[1] for x in s2]
yeast_canda["hitID"] = s20

yeast_canda["combine"] = yeast_canda["geneID"] + "@@" + yeast_canda["hitID"]



from src.mainFunction import *
result_df["pident1"] = singleMapping(canda_yeast["pident"],canda_yeast["combine"],result_df['combine1'])
result_df["pident2"] = singleMapping(yeast_canda["pident"],yeast_canda["combine"],result_df['combine2'])


result_df.to_excel('result/tm_score.xlsx')
g1 = result_df[result_df['tm_score'] < 0.50]
g2 = result_df[result_df['tm_score'] >= 0.50]




# plot for this dataset
result_df['average_pidentity'] = result_df[['pident1', 'pident2']].mean(axis=1)
# drop na
result_df = result_df.dropna(subset=['average_pidentity'])
# density plot
import matplotlib.pyplot as plt
import seaborn as sns
sns.displot(result_df, x="tm_score", height=3, aspect=1.2, alpha=0.25)
plt.xlabel('tm_score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/tm_score_distribution_for_two_yeasts.pdf', bbox_inches='tight')




from scipy import stats
res = stats.pearsonr(result_df['average_pidentity'], result_df['tm_score'])
res



# scatter plot
from scipy.stats import gaussian_kde
from scipy import stats

x = result_df['average_pidentity'].tolist()
y = result_df['tm_score'].tolist()

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
fig, ax = plt.subplots(1,1,figsize=(3.6, 3))
ax.scatter(x, y, c=z, s=20)

plt.xlabel('Pidentity', fontsize=15)
plt.ylabel('tm_score', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
plt.savefig('result/tm_score_&_pidentity_relations_for_two_yeasts.pdf', bbox_inches='tight')





# check the effect of pidentity
result_df_filter1 = result_df[result_df['average_pidentity'] <=30]
result_df_filter1['Pidentity'] = "≤30%"


# calculate the person coefficient
res = stats.pearsonr(result_df_filter1['average_pidentity'], result_df_filter1['tm_score'])
res


result_df_filter2 = result_df[result_df['average_pidentity'] >=70]
result_df_filter2['Pidentity'] = "≥70%"
res = stats.pearsonr(result_df_filter2['average_pidentity'].tolist(), result_df_filter2['tm_score'].tolist())
res


# combine the above two figures
frames = [result_df_filter1, result_df_filter2]
result = pd.concat(frames)
sns.displot(data=result, x='tm_score', hue='Pidentity', fill=True, palette=sns.color_palette('bright')[:2], height=3, aspect=1.2, alpha=0.25,kde=True)
plt.xlabel('tm_score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('result/tm_score_&_pidentity_relations_for_two_yeasts_filter.pdf', bbox_inches='tight')


