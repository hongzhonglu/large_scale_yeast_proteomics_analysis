from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 10

import pandas as pd
import matplotlib.pyplot as plt
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
FIGURE_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'figure')
os.makedirs(FIGURE_DIR, exist_ok=True)
data_path = os.path.join(INTERMEDIATE_DIR, "10fold_R2_results.txt")
df = pd.read_csv(data_path, sep='\t')
n_proteins_manual = 23

folds = sorted(df['fold'].unique())
colors = plt.cm.tab10([(i - 1) / 10 for i in folds])
mean_cv = df.groupby('min_protein_count')['R2_cv'].mean()
mean_test = df.groupby('min_protein_count')['R2_test'].mean()

fig, ax = plt.subplots(figsize=(10, 10), dpi=400)
ax.set_xlabel('Number of Proteins', fontsize=38)
ax.set_ylabel('R²', fontsize=38)
ax.tick_params(axis='y', labelsize=30)
ax.tick_params(axis='x', labelsize=30)
ax.grid(True, linestyle='--', alpha=0.7)
ax.set_xticks(range(0, 65, 5))
ax.set_ylim(0, 1.05)
for i, fold in enumerate(folds):
    sub = df[df['fold'] == fold]
    ax.plot(sub['min_protein_count'], sub['R2_test'], color=colors[i], linestyle='-', label=f'Fold {fold}')
ax.plot(mean_test.index, mean_test.values, color='black', linewidth=2.5, label='Mean')
if n_proteins_manual in mean_test.index:
    r2_test_val = mean_test.loc[n_proteins_manual]
    annotation_text = f"$\\mathbf{{Proteins:}}$ {int(n_proteins_manual)}\n$\\mathbf{{R^2:}}$ {r2_test_val:.3f}"
    annotation_x = n_proteins_manual + 10 if n_proteins_manual <= 30 else n_proteins_manual + 2
    middle_y = r2_test_val * 0.5 + 0.25
    ax.annotate(annotation_text, xy=(n_proteins_manual, r2_test_val), xytext=(annotation_x + 2, middle_y - 0.21),
                 arrowprops=dict(arrowstyle='->', color='black', lw=1.5), fontsize=34,
                 bbox=dict(boxstyle="round,pad=0.4", facecolor='lightyellow', alpha=0.8, edgecolor='gray'),
                 ha='center', va='center')
ax.legend(fontsize=22, loc='lower right', ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "fig4C_R2_test_10fold.png"), dpi=400)
plt.close()
