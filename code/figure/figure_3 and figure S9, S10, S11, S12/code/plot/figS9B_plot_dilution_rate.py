"""
figS9B: Dilution-rate scatter plot (Train=Rosemary, Val=Jianye, Val2=Ibrahim E. Elsemman).
Output: figS9B_plot_dilution_rate_1x1_scatter.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rcParams
import os

rcParams['font.family'] = 'Arial'

# --- Paths (relative) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, 'original_data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)
CSV_PATH = os.path.join(INTERMEDIATE_DIR, 'results_Val_adj_R2_1_8.csv')
PHENOTYPE_PATH = os.path.join(ORIGINAL_DATA_DIR, 'physiology_collection.xlsx')
ORGANELLE_PATH = os.path.join(ORIGINAL_DATA_DIR, 'all_organelle_fraction_test_260412.xlsx')

results_df = pd.read_csv(CSV_PATH)

df_phenotype = pd.read_excel(PHENOTYPE_PATH, sheet_name='Sheet1')
if 'sampleID' in df_phenotype.columns:
    df_phenotype = df_phenotype.set_index('sampleID')
elif df_phenotype.index.name is None and len(df_phenotype.columns) > 0:
    first_col = df_phenotype.columns[0]
    if df_phenotype[first_col].astype(str).str.match(r'^(prot\.|translation|carbon)').any():
        df_phenotype = df_phenotype.set_index(first_col)
df_organelle = pd.read_excel(ORGANELLE_PATH, sheet_name='Sheet1')
df_organelle = df_organelle.iloc[:, 1:]
df_organelle = df_organelle.set_index(df_organelle.columns[0]).T
train_mask = df_phenotype['source'].astype(str).str.contains('Rosemary', case=False, na=False) & df_phenotype.index.astype(str).str.startswith('prot.')
val_mask = df_phenotype['source'].astype(str).str.contains('Jianye', case=False, na=False)
val2_mask = df_phenotype['source'].eq('Ibrahim E. Elsemman')
df_phenotype_train = df_phenotype[train_mask]
df_phenotype_val = df_phenotype[val_mask]
df_phenotype_val2 = df_phenotype[val2_mask]
df_organelle_train = df_organelle.loc[df_organelle.index.intersection(df_phenotype_train.index)]
df_organelle_val = df_organelle.loc[df_organelle.index.intersection(df_phenotype_val.index)]
df_organelle_val2 = df_organelle.loc[df_organelle.index.intersection(df_phenotype_val2.index)]
df_organelle_val2 = df_organelle_val2[~df_organelle_val2.index.astype(str).str.startswith("translation_inhibitor")]
df_phenotype_val2 = df_phenotype_val2[~df_phenotype_val2.index.astype(str).str.startswith("translation_inhibitor")]

target_variable = 'dilution rate (/h)'
exclude_samples = [f'prot.{j}' for j in range(22, 43)]

target_results = results_df[results_df['目标变量'] == target_variable].copy()
param_counts = []
for _, row in target_results.iterrows():
    count = sum(1 for i in range(1, 9) if f'参数{i}' in row and pd.notna(row.get(f'参数{i}')))
    param_counts.append(count)
target_results['特征数量'] = param_counts
candidate = target_results[target_results['特征数量'] == 3]
best_row = candidate.loc[candidate['Val_adj_R2'].idxmax()]
params = []
for k in range(1, 9):
    col = f'参数{k}'
    if col in best_row and pd.notna(best_row.get(col)):
        params.append(str(best_row[col]).replace(' (+)', '').replace(' (-)', ''))
params_raw = list(params)
params = [p for p in params_raw if p in df_organelle_train.columns and p in df_organelle_val.columns and p in df_organelle_val2.columns]
if len(params) < 3:
    params = [p for p in params_raw if p in df_organelle_train.columns and p in df_organelle_val.columns]
params = params[:3]
if 'ribosome' in params:
    params = ['ribosome'] + [p for p in params if p != 'ribosome']


def general_model(X, C, *coefs):
    out = C
    for i in range(X.shape[1]):
        if i < len(coefs):
            out += coefs[i] * X[:, i]
    return out


filtered_organelle_train = df_organelle_train[~df_organelle_train.index.isin(exclude_samples)]
filtered_phenotype_train = df_phenotype_train[~df_phenotype_train.index.isin(exclude_samples)]
common_train = filtered_organelle_train.index.intersection(filtered_phenotype_train.index)
X_train = filtered_organelle_train.loc[common_train][params].values
y_train = filtered_phenotype_train.loc[common_train, target_variable].values
initial_params = [0.0] + [1.0] * len(params)
popt, _ = curve_fit(general_model, X_train, y_train, p0=initial_params)

filtered_organelle_val = df_organelle_val[~df_organelle_val.index.isin(exclude_samples)]
filtered_phenotype_val = df_phenotype_val[~df_phenotype_val.index.isin(exclude_samples)]
common_val = filtered_organelle_val.index.intersection(filtered_phenotype_val.index)
X_val = filtered_organelle_val.loc[common_val][params].values
y_val_true = filtered_phenotype_val.loc[common_val, target_variable].values

common_val2 = df_organelle_val2.index.intersection(df_phenotype_val2.index)
avail_val2 = [c for c in params if c in df_organelle_val2.columns]
X_val2 = df_organelle_val2.loc[common_val2][avail_val2].values if avail_val2 else np.empty((len(common_val2), 0))
y_val2_true = df_phenotype_val2.loc[common_val2, target_variable].values

y_train_pred = general_model(X_train, *popt)
y_val_pred = general_model(X_val, *popt)
if X_val2.shape[1] == len(popt) - 1:
    y_val2_pred = general_model(X_val2, *popt)
else:
    y_val2_pred = np.full(len(y_val2_true), np.nan)

combined_color_train = '#FF6347'
combined_color_val = '#FFFF00'
color_ibrahim = '#87CEEB'
plt.figure(figsize=(32, 8), dpi=400)
plt.scatter(y_train, y_train_pred, color=combined_color_train, edgecolor='darkblue', alpha=0.7, s=120, label='Rosemary')
plt.scatter(y_val_true, y_val_pred, color=combined_color_val, edgecolor='darkblue', alpha=0.7, s=120, label='Jianye')
if not np.all(np.isnan(y_val2_pred)):
    plt.scatter(y_val2_true, y_val2_pred, color=color_ibrahim, edgecolor='darkblue', alpha=0.7, s=120, label='Ibrahim E. Elsemman')
all_x = np.concatenate([y_train, y_val_true, y_val2_true])
all_y = np.concatenate([y_train_pred, y_val_pred, y_val2_pred if not np.all(np.isnan(y_val2_pred)) else y_val2_true])
mn = min(np.min(all_x), np.min(all_y))
mx = max(np.max(all_x), np.max(all_y))
line_vals = np.linspace(mn, mx, 100)
plt.plot(line_vals, line_vals, 'k--', linewidth=1, alpha=0.7)
plt.xlabel('Measured value\n(/h)', fontsize=34)
plt.ylabel('Predicted value\n(/h)', fontsize=34)
plt.title('Dilution rate', fontsize=38, fontweight='bold', pad=15)
plt.tick_params(axis='both', labelsize=28)
plt.legend(loc='lower right', fontsize=28)
plt.xlim([mn, mx])
plt.ylim([mn, mx])
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
n_train = len(y_train)
n_val = len(y_val_true)
p_num = len(params)
r2_train = r2_score(y_train, y_train_pred)
r2_val = r2_score(y_val_true, y_val_pred)
adj_r2_train = 1 - (1 - r2_train) * (n_train - 1) / max(n_train - p_num - 1, 1)
adj_r2_val = 1 - (1 - r2_val) * (n_val - 1) / max(n_val - p_num - 1, 1)
if not np.all(np.isnan(y_val2_pred)) and len(avail_val2) > 0:
    n_val2 = len(y_val2_true)
    p_val2 = len(avail_val2)
    r2_val2 = r2_score(y_val2_true, y_val2_pred)
    adj_r2_val2 = 1 - (1 - r2_val2) * (n_val2 - 1) / max(n_val2 - p_val2 - 1, 1)
else:
    adj_r2_val2 = np.nan
line_height = 0.085
y_pos = 0.95
plt.text(0.05, y_pos, r"$\bf{C}$" + f" = {popt[0]:.3f}", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
y_pos -= line_height
greek_letters = ['α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ']
for j in range(min(len(params), len(popt) - 1)):
    pname = params[j] if j < len(params) else f'p{j}'
    pval = popt[j + 1]
    if len(pname) > 25:
        plt.text(0.05, y_pos, r"$\bf{" + greek_letters[j] + r"}$" + f" = {pval:.3f}", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
        y_pos -= line_height
        plt.text(0.05, y_pos, f"({pname})", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
    else:
        plt.text(0.05, y_pos, r"$\bf{" + greek_letters[j] + r"}$" + f" = {pval:.3f} ({pname})", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
    y_pos -= line_height
y_pos -= 0.01
plt.text(0.05, y_pos, r"$\bf{adj.R^2}\ (\bf{R})$ = " + f"{adj_r2_train:.3f}", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
y_pos -= line_height
plt.text(0.05, y_pos, r"$\bf{adj.R^2}\ (\bf{J})$ = " + f"{adj_r2_val:.3f}", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')
y_pos -= line_height
plt.text(0.05, y_pos, r"$\bf{adj.R^2}\ (\bf{I})$ = " + f"{adj_r2_val2:.3f}", transform=plt.gca().transAxes, fontsize=28, verticalalignment='top')

plt.tight_layout()
out_path = os.path.join(FIGURE_DIR, 'figS9B_plot_dilution_rate_1x1_scatter.png')
plt.savefig(out_path, dpi=400, bbox_inches='tight')
plt.show()
