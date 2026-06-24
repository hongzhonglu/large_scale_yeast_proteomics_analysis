import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
OUTPUT_XLSX_PATH = os.path.join(INTERMEDIATE_DIR, 'all_summary.xlsx')
OUTPUT_CSV_PATH = os.path.join(INTERMEDIATE_DIR, 'results_Val_adj_R2_1_8.csv')


def add_sign(param_series, param_value):
    values = []
    for i, param in enumerate(param_series):
        if pd.isna(param):
            values.append(param)
            continue
        pname = str(param).replace(' (+)', '').replace(' (-)', '')
        v = param_value.iloc[i]
        sign = '+' if pd.notnull(v) and v > 0 else '-' if pd.notnull(v) and v < 0 else ''
        values.append(f'{pname} ({sign})')
    return pd.Series(values, index=param_series.index)


def main():
    # Read the existing merged results CSV and only enrich with (+)/(-) signs.
    # (Per-target CSVs are no longer produced; the fit script writes one combined CSV.)
    if not os.path.exists(OUTPUT_CSV_PATH):
        print(f'not found: {OUTPUT_CSV_PATH}')
        return
    results_df = pd.read_csv(OUTPUT_CSV_PATH)
    if 'Val_adj_R2' in results_df.columns:
        results_df = results_df.sort_values('Val_adj_R2', ascending=False).reset_index(drop=True)
    for i in range(1, 21):
        param_col = f'参数{i}'
        param_value_col = f'param{i}'
        if param_col in results_df.columns and param_value_col in results_df.columns:
            results_df[param_col] = add_sign(results_df[param_col], results_df[param_value_col])
    results_df.to_csv(OUTPUT_CSV_PATH, index=False)
    results_df.to_excel(OUTPUT_XLSX_PATH, index=False)
    print(f'saved: {OUTPUT_CSV_PATH}')
    print(f'saved: {OUTPUT_XLSX_PATH}')


if __name__ == '__main__':
    main()
