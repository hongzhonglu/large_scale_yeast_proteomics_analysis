"""
fig5A & fig5B: Best-formula scatter plots (Train and Val) for the 4 target variables.
- fig5A: plot_best_scatter_train_1x4.png  (Train set, Rosemary)
- fig5B: plot_best_scatter_val_1x4.png    (Val set, Jianye)
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr, spearmanr
from matplotlib import rcParams
from scipy.optimize import curve_fit
import os

rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 12

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
df_phenotype_train = df_phenotype[train_mask]
df_phenotype_val = df_phenotype[val_mask]
df_organelle_train = df_organelle.loc[df_organelle.index.intersection(df_phenotype_train.index)]
df_organelle_val = df_organelle.loc[df_organelle.index.intersection(df_phenotype_val.index)]

target_variables = ['dilution rate (/h)', 'qGlucose (mmol/gDW h)', 'qCO2 (mmol/gDW h)', 'qO2 (mmol/gDW h)']

inflection_points = {
    'dilution rate (/h)': 3,
    'qGlucose (mmol/gDW h)': 2,
    'qCO2 (mmol/gDW h)': 2,
    'qO2 (mmol/gDW h)': 2
}


def general_model(X, C, *params):
    result = C
    for i in range(X.shape[1]):
        if i < len(params):
            result += params[i] * X[:, i]
    return result


best_models = []
for target in target_variables:
    best_feature_count = inflection_points[target]
    target_results = results_df[results_df['目标变量'] == target].copy()
    param_counts = []
    for _, row in target_results.iterrows():
        count = 0
        for i in range(1, 9):
            col = f'参数{i}'
            if col in row and pd.notna(row[col]):
                count += 1
        param_counts.append(count)
    target_results['特征数量'] = param_counts
    candidate = target_results[target_results['特征数量'] == best_feature_count]
    if candidate.empty:
        continue
    best_row = candidate.loc[candidate['Val_adj_R2'].idxmax()]
    best_models.append(best_row)

best_models = pd.DataFrame(best_models).reset_index(drop=True)

combined_color_train = '#FF6347'
combined_color_val = '#FFFF00'
formula_on_x_axis = False
train_panels = []
val_panels = []
train_axis_limits_by_target = {}

for i, (_, model) in enumerate(best_models.iterrows()):
    if i >= 4:
        break
    target_variable = model['目标变量']
    params = []
    for k in range(1, 9):
        col = f'参数{k}'
        if col in model and pd.notna(model[col]):
            name = str(model[col]).replace(' (+)', '').replace(' (-)', '')
            params.append(name)
    params = [p for p in params if p in df_organelle_train.columns and p in df_organelle_val.columns]
    if 'ribosome' in params:
        params = ['ribosome'] + [p for p in params if p != 'ribosome']
    num_features = len(params)
    if num_features == 0:
        continue
    param_values = []
    combined_c = 0.0
    if target_variable == 'dilution rate (/h)':
        exclude_samples = [f'prot.{j}' for j in range(22, 43)]
        filtered_train_organelle = df_organelle_train[~df_organelle_train.index.isin(exclude_samples)]
        filtered_train_phenotype = df_phenotype_train[~df_phenotype_train.index.isin(exclude_samples)]
        common_train = filtered_train_organelle.index.intersection(filtered_train_phenotype.index)
        X_train = filtered_train_organelle.loc[common_train][params].values
        y_true_train = filtered_train_phenotype.loc[common_train, target_variable].values
    else:
        common_train = df_organelle_train.index.intersection(df_phenotype_train.index)
        X_train = df_organelle_train.loc[common_train][params].values
        y_true_train = df_phenotype_train.loc[common_train, target_variable].values
    initial_params = [0.0] + [1.0] * num_features
    popt, pcov = curve_fit(general_model, X_train, y_true_train, p0=initial_params)
    combined_c = popt[0]
    param_values = list(popt[1:])
    y_pred_train = general_model(X_train, *popt)
    r2_train = r2_score(y_true_train, y_pred_train)
    n_train = len(y_true_train)
    p_train = num_features
    train_adj_r2 = 1 - (1 - r2_train) * (n_train - 1) / (n_train - p_train - 1)
    train_rmse = np.sqrt(mean_squared_error(y_true_train, y_pred_train))
    try:
        train_pcc, train_pcc_p = pearsonr(y_true_train, y_pred_train)
    except Exception:
        train_pcc, train_pcc_p = np.nan, np.nan
    try:
        train_scc, train_scc_p = spearmanr(y_true_train, y_pred_train)
    except Exception:
        train_scc, train_scc_p = np.nan, np.nan
    if target_variable == 'dilution rate (/h)':
        exclude_samples = [f'prot.{j}' for j in range(22, 43)]
        filtered_val_organelle = df_organelle_val[~df_organelle_val.index.isin(exclude_samples)]
        filtered_val_phenotype = df_phenotype_val[~df_phenotype_val.index.isin(exclude_samples)]
        common_val = filtered_val_organelle.index.intersection(filtered_val_phenotype.index)
        X_val = filtered_val_organelle.loc[common_val][params].values
        y_true_val = filtered_val_phenotype.loc[common_val, target_variable].values
    else:
        common_val = df_organelle_val.index.intersection(df_phenotype_val.index)
        X_val = df_organelle_val.loc[common_val][params].values
        y_true_val = df_phenotype_val.loc[common_val, target_variable].values
    if target_variable == 'qGlucose (mmol/gDW h)':
        y_true_val = np.abs(y_true_val)
    y_pred_val = general_model(X_val, *popt)
    r2_val = r2_score(y_true_val, y_pred_val)
    n_val = len(y_true_val)
    p_val = num_features
    val_adj_r2 = 1 - (1 - r2_val) * (n_val - 1) / (n_val - p_val - 1)
    val_rmse = np.sqrt(mean_squared_error(y_true_val, y_pred_val))
    try:
        val_pcc, val_pcc_p = pearsonr(y_true_val, y_pred_val)
    except Exception:
        val_pcc, val_pcc_p = np.nan, np.nan
    try:
        val_scc, val_scc_p = spearmanr(y_true_val, y_pred_val)
    except Exception:
        val_scc, val_scc_p = np.nan, np.nan

    def format_axis_label(base_label, target):
        if target == 'dilution rate (/h)':
            return f"{base_label}\n(/h)"
        elif target == 'qGlucose (mmol/gDW h)':
            return f"{base_label}\n(mmol/gDW h)"
        elif target == 'qCO2 (mmol/gDW h)':
            return f"{base_label}\n(mmol/gDW h)"
        elif target == 'qO2 (mmol/gDW h)':
            return f"{base_label}\n(mmol/gDW h)"
        return base_label

    def get_title_without_unit(target):
        if target == 'dilution rate (/h)':
            return 'Dilution rate'
        elif target == 'qGlucose (mmol/gDW h)':
            return 'qGlucose'
        elif target == 'qCO2 (mmol/gDW h)':
            return 'qCO2'
        elif target == 'qO2 (mmol/gDW h)':
            return 'qO2'
        return target

    if formula_on_x_axis:
        plot_x_train, plot_y_train = y_pred_train, y_true_train
        xlabel = format_axis_label("Predicted value", target_variable)
        ylabel = format_axis_label("Measured value", target_variable)
    else:
        plot_x_train, plot_y_train = y_true_train, y_pred_train
        xlabel = format_axis_label("Measured value", target_variable)
        ylabel = format_axis_label("Predicted value", target_variable)
    train_panels.append({
        'target_variable': get_title_without_unit(target_variable),
        'params': params,
        'param_values': param_values,
        'combined_c': combined_c,
        'x': plot_x_train,
        'y': plot_y_train,
        'xlabel': xlabel,
        'ylabel': ylabel,
        'adj_r2': train_adj_r2,
        'rmse': train_rmse,
        'pcc': train_pcc,
        'scc': train_scc
    })

    if formula_on_x_axis:
        plot_x_val, plot_y_val = y_pred_val, y_true_val
    else:
        plot_x_val, plot_y_val = y_true_val, y_pred_val
    val_panels.append({
        'target_variable': get_title_without_unit(target_variable),
        'params': params,
        'param_values': param_values,
        'combined_c': combined_c,
        'x': plot_x_val,
        'y': plot_y_val,
        'xlabel': xlabel,
        'ylabel': ylabel,
        'adj_r2': val_adj_r2,
        'rmse': val_rmse,
        'pcc': val_pcc,
        'scc': val_scc
    })

# --- fig5A: Train scatter ---
fig1, axes1 = plt.subplots(1, 4, figsize=(16, 5), dpi=400)
axes1 = axes1.flatten()
for i, panel in enumerate(train_panels[:4]):
    ax = axes1[i]
    ax.scatter(panel['x'], panel['y'], color=combined_color_train, edgecolor='darkblue', alpha=0.7, s=20, label='Data Points')
    min_val = min(np.min(panel['x']), np.min(panel['y']))
    max_val = max(np.max(panel['x']), np.max(panel['y']))
    line_vals = np.linspace(min_val, max_val, 100)
    ax.plot(line_vals, line_vals, 'k--', linewidth=1, alpha=0.7)
    ax.set_title(f"{panel['target_variable']}", fontsize=18, fontweight='bold')
    y_pos = 0.95
    line_height = 0.05
    ax.text(0.05, y_pos, f"C = {panel['combined_c']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    greek_letters = ['α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ']
    for j, (param_name, param_value) in enumerate(zip(panel['params'], panel['param_values'])):
        if j < len(greek_letters):
            if len(param_name) > 25:
                ax.text(0.05, y_pos, f"{greek_letters[j]} = {param_value:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
                y_pos -= line_height
                ax.text(0.05, y_pos, f"({param_name})", transform=ax.transAxes, fontsize=10, verticalalignment='top')
            else:
                ax.text(0.05, y_pos, f"{greek_letters[j]} = {param_value:.3f} ({param_name})", transform=ax.transAxes, fontsize=10, verticalalignment='top')
            y_pos -= line_height
    ax.text(0.05, y_pos, f"\nadj.R² = {panel['adj_r2']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nRMSE = {panel['rmse']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nPCC = {panel['pcc']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nSCC = {panel['scc']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    ax.set_xlabel(panel['xlabel'], fontsize=20)
    ax.set_ylabel(panel['ylabel'], fontsize=20)
    if i == 0:
        ax.legend(fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=16)
    axis_limits = [min_val, max_val]
    train_axis_limits_by_target[panel['target_variable']] = axis_limits
    ax.set_xlim(axis_limits)
    ax.set_ylim(axis_limits)
    ax.set_aspect('equal')

plt.tight_layout()
out_png1 = os.path.join(FIGURE_DIR, 'fig4A_plot_best_scatter_train_1x4.png')
plt.savefig(out_png1, dpi=400, bbox_inches='tight')
plt.show()

# --- fig5B: Val scatter ---
fig2, axes2 = plt.subplots(1, 4, figsize=(16, 5), dpi=400)
axes2 = axes2.flatten()
for i, panel in enumerate(val_panels[:4]):
    ax = axes2[i]
    ax.scatter(panel['x'], panel['y'], color=combined_color_val, edgecolor='darkblue', alpha=0.7, s=20, label='Data Points')
    min_val = min(np.min(panel['x']), np.min(panel['y']))
    max_val = max(np.max(panel['x']), np.max(panel['y']))
    line_vals = np.linspace(min_val, max_val, 100)
    ax.plot(line_vals, line_vals, 'k--', linewidth=1, alpha=0.7)
    ax.set_title(f"{panel['target_variable']}", fontsize=18, fontweight='bold')
    y_pos = 0.95
    line_height = 0.05
    ax.text(0.05, y_pos, f"C = {panel['combined_c']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    greek_letters = ['α', 'β', 'γ', 'δ', 'ε', 'ζ', 'η', 'θ', 'ι', 'κ']
    for j, (param_name, param_value) in enumerate(zip(panel['params'], panel['param_values'])):
        if j < len(greek_letters):
            if len(param_name) > 25:
                ax.text(0.05, y_pos, f"{greek_letters[j]} = {param_value:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
                y_pos -= line_height
                ax.text(0.05, y_pos, f"({param_name})", transform=ax.transAxes, fontsize=10, verticalalignment='top')
            else:
                ax.text(0.05, y_pos, f"{greek_letters[j]} = {param_value:.3f} ({param_name})", transform=ax.transAxes, fontsize=10, verticalalignment='top')
            y_pos -= line_height
    ax.text(0.05, y_pos, f"\nadj.R² = {panel['adj_r2']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nRMSE = {panel['rmse']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nPCC = {panel['pcc']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    y_pos -= line_height
    ax.text(0.05, y_pos, f"\nSCC = {panel['scc']:.3f}", transform=ax.transAxes, fontsize=10, verticalalignment='top')
    ax.set_xlabel(panel['xlabel'], fontsize=20)
    ax.set_ylabel(panel['ylabel'], fontsize=20)
    if i == 0:
        ax.legend(fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=16)
    axis_limits = train_axis_limits_by_target.get(panel['target_variable'], [min_val, max_val])
    ax.set_xlim(axis_limits)
    ax.set_ylim(axis_limits)
    ax.set_aspect('equal')

plt.tight_layout()
out_png2 = os.path.join(FIGURE_DIR, 'fig4B_plot_best_scatter_val_1x4.png')
plt.savefig(out_png2, dpi=400, bbox_inches='tight')
plt.show()
