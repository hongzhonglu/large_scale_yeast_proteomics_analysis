# orthlogFinder process
import pandas as pd
import os


def getAlltmScore(data_dir, group):
    all_file = os.listdir(data_dir)
    pair = []
    tm_score = []
    identity_us_align = []
    # process:
    for x in all_file:
        print(x)
        file01 = open(us_align_dir0 + x).readlines()
        if len(file01) >= 2:
            tm_score1 = file01[14]
            tm_score2 = file01[15]
            identity = float(file01[13].split("= ")[3])
            v1 = float(tm_score1.split(' (')[0].split('= ')[1])
            v2 = float(tm_score2.split(' (')[0].split('= ')[1])
            mean0 = (v1 + v2) / 2
            pair.append(x)
            tm_score.append(mean0)
            identity_us_align.append(identity)
    def find_between_r(s, first, last):
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
    result_df = pd.DataFrame({"pair_id": pair, "tm_score": tm_score, "identity_us_align": identity_us_align})
    result_df1 = result_df['pair_id'].str.split('@@', 1, expand=True)
    result_df1.columns = ['pro1', 'pro2']
    result_df2 = result_df1.copy()
    result_df2['pro1'] = find_between_r(s=result_df1['pro1'].tolist(), first='AF-', last='-F1')
    result_df2['pro2'] = find_between_r(s=result_df1['pro2'].tolist(), first='AF-', last='-F1')
    result_df['pro1'] = result_df2['pro1']
    result_df['pro2'] = result_df2['pro2']
    result_df['group'] = group
    return result_df



## part1
group_name = 'emp_vs_emp'
us_align_dir0 = '/Users/xluhon/Documents/data_for_structure_align/tm_score_subpathway/' + group_name + '/'
df_emp_vs_emp = getAlltmScore(data_dir=us_align_dir0, group=group_name)


group_name = 'emp_vs_tca'
us_align_dir0 = '/Users/xluhon/Documents/data_for_structure_align/tm_score_subpathway/' + group_name + '/'
df_emp_vs_tca = getAlltmScore(data_dir=us_align_dir0, group=group_name)

group_name = 'emp_vs_ppp'
us_align_dir0 = '/Users/xluhon/Documents/data_for_structure_align/tm_score_subpathway/' + group_name + '/'
df_emp_vs_ppp = getAlltmScore(data_dir=us_align_dir0, group=group_name)

group_name = 'emp_vs_op'
us_align_dir0 = '/Users/xluhon/Documents/data_for_structure_align/tm_score_subpathway/' + group_name + '/'
df_emp_vs_op = getAlltmScore(data_dir=us_align_dir0, group=group_name)


# merge two dataframe
frames = [df_emp_vs_emp, df_emp_vs_tca, df_emp_vs_ppp, df_emp_vs_op]
result = pd.concat(frames)

result = result[result['tm_score'] <= 0.6]

# boxplot
import matplotlib.pyplot as plt
import seaborn as sns
sns.boxplot(data=result, x="group", y="tm_score")








## part2
"""
# this script is running on the linux as it takes a lot of time in calculation.
group_name = 'pro_with_unknown_function'
us_align_dir0 = '/home/yeast/Documents/tm_out/'
df = getAlltmScore(data_dir=us_align_dir0, group=group_name)
df.to_csv("/home/yeast/Documents/tm_score_for_pro_with_unknown_function.csv")
"""


# filter
df_all = pd.read_csv("/Users/xluhon/Documents/data_for_structure_align/tm_score_for_pro_with_unknown_function.csv")
df_all_filter = df_all[df_all["tm_score"] >=0.5]
df_all_filter.to_csv("result/tm_score_for_pro_with_unknown_function.csv")

df_all_filter_other = df_all[df_all["tm_score"] < 0.5]
# compare proteins with tm_score larger or lower than 0.5
n1 = len(set(df_all_filter['pro1'].tolist()))
n2 = len(set(df_all_filter_other['pro1'].tolist())-set(df_all_filter['pro1'].tolist()))



# plot
sns.displot(df_all_filter, x="tm_score")
plt.xlabel('tm_score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

sns.displot(df_all, x="tm_score")
plt.xlabel('tm_score', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


sns.displot(df_all_filter, x="identity_us_align")
plt.xlabel('identity_us_align', fontsize=15)
plt.ylabel('Count', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)



# scatter plot
from scipy.stats import gaussian_kde
import numpy as np
from scipy import stats

x = df_all_filter['identity_us_align'].tolist()
y = df_all_filter['tm_score'].tolist()

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=20)

plt.xlabel('identity_us_align', fontsize=15)
plt.ylabel('tm_score', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()
res = stats.pearsonr(df_all_filter['identity_us_align'].tolist(), df_all_filter['tm_score'].tolist())
res



# in-depth analysis of these proteins with tm-score larger than 0.5
# pro1: uncharacterized function
# pro2: characterized function
# Uncharacterized protein, Putative uncharacterized protein
# Q04516 · AIM33_YEAST
df_all_filter = df_all_filter.sort_values('tm_score', ascending=False)



