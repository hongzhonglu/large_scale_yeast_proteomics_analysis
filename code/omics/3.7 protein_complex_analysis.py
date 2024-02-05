# Note:
# all these analysis is based on the protein abundance in the unit of mmol/DCW


# import self function
from src.model_process import *
from src.protein_process import *
import matplotlib.pyplot as plt
import seaborn as sns

def find_outlier(df, column):
    # Find first and third quartile
    q1 = df[column].quantile(0.25)
    q3 = df[column].quantile(0.75)

    # Find interquartile range
    IQR = q3 - q1

    # Find lower and upper bound
    lower_bound = q1 - 1.5 * IQR
    upper_bound = q3 + 1.5 * IQR

    # Remove outliers, here actuall the too larger or too smaller was replaced by NAN.
    df[column] = df[column][df[column] > lower_bound]
    df[column] = df[column][df[column] < upper_bound]

    return df[column]


# data preprocess
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx") # the unit the g/gDW???? should be wrong


# analyze the protein resouce allocation at stoichiometric level
protein_complex = pd.read_excel("data/complex_info.xlsx")
protein_complex["subunit"] = protein_complex["subunit"].str.replace("-MONOMER","")

# combine this with the protein abundance information
result = pd.merge(protein_complex, omics_combine_auto, left_on='subunit', right_on='all_gene', how='left')
result.to_excel("data/proteomics/protein_complex_subunit_with_abundance.xlsx")

# use a function to run the above task
def complex_subunit_analysis(complexid, min_num=10):
    test_data = result[result['complex'] == complexid]
    test_data0 = test_data[["subunit"]]
    test_data1 = test_data.iloc[:, 5:]
    test_data1 = test_data1.iloc[:, 0:].div(test_data.iloc[:, 3], axis=0) # remove the effect of coefficient
    col0 = list(test_data1.columns)
    # loop for the column
    for i in col0:
        print(i)
        ss = test_data1[i].tolist()
        len_over_zero = len([x for x in ss if x > 0])
        if len_over_zero > 0.5 * len(ss):
            test_data0[i] = ss
    test_data_x = test_data0.copy().iloc[:, 1:]
    test_data_x['Total'] = test_data_x.count(axis=1)
    gene_id = test_data0['subunit'].tolist()
    gene_count = test_data_x['Total'].tolist()
    gene_id_select = [x for x, y in zip(gene_id, gene_count) if y >= 10]
    test_data0 = test_data0[test_data0['subunit'].isin(gene_id_select)]
    if test_data0.shape[1] >= min_num and test_data0.shape[0] >= 2:
        # graph plot step
        # bar plot
        title0 = 'result/figure/' + complexid + '.pdf'
        print(title0)
        plt.figure(figsize=(4, 4))
        sns.barplot(x="subunit", y="count", data=test_data, capsize=.2)
        plt.xlabel(complexid + "_subunit", fontsize=12)
        plt.ylabel("Count", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(rotation=45)
        plt.savefig(title0, bbox_inches='tight')
        # box plot
        col1 = list(test_data0.columns)
        # loop for the column
        s1 = []
        s2 = []
        for i in col1[2:]:
            print(i)
            s1 = s1 + test_data0['subunit'].tolist()
            s2 = s2 + test_data0[i].tolist()
        new_df = pd.DataFrame({"Subunit": s1, "Abundance": s2})

        new_df['Abundance'] = find_outlier(new_df, 'Abundance') # remove outlier data point

        title1 = 'result/figure/' + complexid + '_subunit_abundance.pdf'

        plt.figure(figsize=(4, 4))
        sns.boxplot(data=new_df, x="Subunit", y="Abundance", hue="Subunit", dodge=False)
        plt.xlabel(complexid + "_subunit", fontsize=12)
        plt.ylabel("Count", fontsize=15)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.xticks(rotation=45)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        plt.savefig(title1, bbox_inches='tight')
        # carry out statistical analysis
        from itertools import combinations
        from scipy.stats import ttest_ind

        new_df = new_df.dropna()
        grps = new_df['Subunit'].unique()
        combs = combinations(grps, 2)
        ttests = {
            f'{c1}_{c2}': ttest_ind(
            new_df.loc[new_df['Subunit'] == c1, 'Abundance'],
            new_df.loc[new_df['Subunit'] == c2, 'Abundance']
            ) for c1, c2 in combs
            }
        print(ttests)
        pvalue = []
        for x, y in ttests.items():
            print(x, y)
            pvalue0 = y[1]
            pvalue.append(pvalue0)
        if min(pvalue) > 0.05:
            return "no_difference"
        else:
            return "with_difference"
    else:
        return "no_valuable_data"


# test the function
complexid ="CPX-1153" #"CPX-568"
complexid ="CPX-568"
complexid ="CPX-1898"
complexid ="CPX-1631"

CPX-1674
CPX-3215
CPX-1165
complexid ="CPX-3088"
ss = complex_subunit_analysis(complexid)


all_complex_id = list(set(protein_complex['complex'].tolist()))
out = []
for xx in all_complex_id:
    print(xx)
    ss0 = complex_subunit_analysis(complexid=xx)
    out.append(ss0)

combine_df = pd.DataFrame({"complexid":all_complex_id,"statistical_analysis":out})
combine_df.to_excel("data/proteomics/protein_complex_subunit_abundance_statistical_analysis.xlsx")











