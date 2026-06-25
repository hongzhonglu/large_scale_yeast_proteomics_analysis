"""
figS9A: Train/Val adjusted R^2 vs number of features (organelles), 1x4 subplots.
Output: figS9A_plot_Train_Val_adjR2_1x4.png
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rcParams

rcParams['font.family'] = 'Arial'

# --- Paths (relative) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)
CSV_PATH = os.path.join(INTERMEDIATE_DIR, 'results_Val_adj_R2_1_8.csv')

results_df = pd.read_csv(CSV_PATH)

target_variables = ['dilution rate (/h)', 'qGlucose (mmol/gDW h)', 'qCO2 (mmol/gDW h)', 'qO2 (mmol/gDW h)']

plot_data = {}
for target in target_variables:
    target_df = results_df[results_df['目标变量'] == target].copy()
    target_df['特征数量'] = target_df.apply(lambda row: sum(pd.notnull(row.get(f'参数{i}', np.nan)) for i in range(1, 9)), axis=1)
    if target not in plot_data:
        plot_data[target] = {
            '特征数量': [],
            '最佳Train_adj_R2': [],
            '最佳Val_adj_R2': [],
            '最佳参数组合': [],
            '拐点用Train_adj_R2': []
        }
    for feature_count in range(1, 9):
        subset = target_df[target_df['特征数量'] == feature_count]
        if not subset.empty:
            best_train_row = subset.loc[subset['Train_adj_R2'].idxmax()]
            best_val_row = subset.loc[subset['Val_adj_R2'].idxmax()]
            plot_data[target]['特征数量'].append(feature_count)
            plot_data[target]['最佳Train_adj_R2'].append(best_train_row['Train_adj_R2'])
            plot_data[target]['最佳Val_adj_R2'].append(best_val_row['Val_adj_R2'])
            plot_data[target]['拐点用Train_adj_R2'].append(best_val_row['Train_adj_R2'])
            params = []
            for j in range(1, 9):
                param_col = f'参数{j}'
                if param_col in best_val_row and pd.notnull(best_val_row.get(param_col)):
                    params.append(str(best_val_row[param_col]))
            plot_data[target]['最佳参数组合'].append(params)

inflection_points = {
    'dilution rate (/h)': 3,
    'qGlucose (mmol/gDW h)': 2,
    'qCO2 (mmol/gDW h)': 2,
    'qO2 (mmol/gDW h)': 2
}

fig, axes = plt.subplots(1, 4, figsize=(32, 8), dpi=400)
axes = axes.flatten()

for i, target in enumerate(target_variables):
    ax = axes[i]
    ax.set_xlabel('Organelles count', fontsize=34)
    if i == 0:
        ax.set_ylabel('Adj.r²', fontsize=34)
    ax.plot(plot_data[target]['特征数量'], plot_data[target]['最佳Train_adj_R2'],
            marker='o', linestyle='-', color='tab:blue', label='Train. adj.R²')
    ax.plot(plot_data[target]['特征数量'], plot_data[target]['最佳Val_adj_R2'],
            marker='s', linestyle='--', color='tab:red', label='Val. adj.R²')
    ax.tick_params(axis='y', labelsize=26)
    ax.tick_params(axis='x', labelsize=26)
    ax.set_ylim([0, 1.05])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xticks(range(1, 9))
    ax.set_xlim([0.5, 8.5])

    def get_title_without_unit(target):
        if target == 'dilution rate (/h)':
            return 'Dilution rate\n(/h)'
        elif target == 'qGlucose (mmol/gDW h)':
            return 'qGlucose\n(mmol/gDW h)'
        elif target == 'qCO2 (mmol/gDW h)':
            return 'qCO2\n(mmol/gDW h)'
        elif target == 'qO2 (mmol/gDW h)':
            return 'qO2\n(mmol/gDW h)'
        return target

    ax.set_title(get_title_without_unit(target), fontsize=34, fontweight='bold', pad=15)
    if i == 0:
        ax.legend(loc='lower left', fontsize=22)

    if target in inflection_points:
        inflection_point = inflection_points[target]
        if inflection_point in plot_data[target]['特征数量']:
            idx = plot_data[target]['特征数量'].index(inflection_point)
            train_r2 = plot_data[target]['拐点用Train_adj_R2'][idx]
            val_r2 = plot_data[target]['最佳Val_adj_R2'][idx]
            params = plot_data[target]['最佳参数组合'][idx]
            ribosome_first = [p for p in params if str(p).replace(' (+)', '').replace(' (-)', '') == 'ribosome']
            others = [p for p in params if str(p).replace(' (+)', '').replace(' (-)', '') != 'ribosome']
            params = ribosome_first + others
            param_str = '\n'.join(params)
            annotation_text = (f"$\\mathbf{{Parameters:}}$ \n{param_str}\n"
                               f"$\\mathbf{{Val.\\ R^2:}}$ {val_r2:.3f}\n"
                               f"$\\mathbf{{Train.\\ R^2:}}$ {train_r2:.3f}")

            if target == 'dilution rate (/h)':
                middle_y = 0.55
                annotation_x = inflection_point - 0.3
            elif target == 'qGlucose (mmol/gDW h)':
                middle_y = 0.3
                annotation_x = inflection_point + 1.3
            elif target == 'qCO2 (mmol/gDW h)':
                middle_y = 0.3
                annotation_x = inflection_point + 0.8
            elif target == 'qO2 (mmol/gDW h)':
                middle_y = 0.3
                annotation_x = inflection_point - 0.4
            else:
                middle_y = 0.3
                annotation_x = inflection_point + 2

            ax.annotate(annotation_text,
                        xy=(inflection_point, val_r2),
                        xytext=(annotation_x, middle_y),
                        arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
                        fontsize=24,
                        bbox=dict(boxstyle="round,pad=0.4", facecolor='white', alpha=0.8, edgecolor='gray'),
                        ha='left', va='center')

for i in range(len(target_variables), len(axes)):
    axes[i].axis('off')

plt.tight_layout()
out_path = os.path.join(FIGURE_DIR, 'figS9A_plot_Train_Val_adjR2_1x4.png')
plt.savefig(out_path, dpi=400, bbox_inches='tight')
plt.show()
