import os
import math
import heapq
import multiprocessing
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from tqdm import tqdm

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, 'original_data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
ORGANELLE_ALL_PATH = os.path.join(ORIGINAL_DATA_DIR, 'main_organelle_fraction_260414.xlsx')
PHYSIOLOGY_PATH = os.path.join(ORIGINAL_DATA_DIR, 'physiology_collection.xlsx')
#旧版是physiology_collection.xlsx
OUTPUT_CSV_PATH = os.path.join(INTERMEDIATE_DIR, 'results_Val_adj_R2_1_8.csv')
OUTPUT_XLSX_PATH = os.path.join(INTERMEDIATE_DIR, 'all_summary.xlsx')
TARGET_VARIABLES = ['dilution rate (/h)', 'qGlucose (mmol/gDW h)', 'qCO2 (mmol/gDW h)', 'qO2 (mmol/gDW h)']
MIN_FEATURES = 1
MAX_FEATURES = 8
TOP_K_PER_FEATURE_COUNT = 10
MAX_WORKERS = min(max(multiprocessing.cpu_count() - 1, 1), 24)

_X_TRAIN = None
_Y_TRAIN = None
_X_VAL = None
_Y_VAL = None
_FEATURE_NAMES = None


def general_model(X, C, *params):
    out = C
    for i in range(X.shape[1]):
        if i < len(params):
            out += params[i] * X[:, i]
    return out


def adjusted_r2(y_true, y_pred, p):
    n = len(y_true)
    if n <= p + 1:
        return np.nan
    r2 = r2_score(y_true, y_pred)
    return 1 - (1 - r2) * (n - 1) / (n - p - 1)


def _init_worker(x_train, y_train, x_val, y_val, feature_names):
    global _X_TRAIN, _Y_TRAIN, _X_VAL, _Y_VAL, _FEATURE_NAMES
    _X_TRAIN = x_train
    _Y_TRAIN = y_train
    _X_VAL = x_val
    _Y_VAL = y_val
    _FEATURE_NAMES = feature_names


def _fit_combo(combo):
    try:
        idx = list(combo)
        x_train = _X_TRAIN[:, idx]
        x_val = _X_VAL[:, idx]
        p0 = [0.0] + [1.0] * len(idx)
        popt, _ = curve_fit(general_model, x_train, _Y_TRAIN, p0=p0, maxfev=50000)
        y_train_pred = general_model(x_train, *popt)
        y_val_pred = general_model(x_val, *popt)
        train_adj = adjusted_r2(_Y_TRAIN, y_train_pred, len(idx))
        val_adj = adjusted_r2(_Y_VAL, y_val_pred, len(idx))
        result = {
            'Train_adj_R2': train_adj,
            'Val_adj_R2': val_adj,
            'C': float(popt[0]),
        }
        for i, cidx in enumerate(idx, 1):
            result[f'参数{i}'] = _FEATURE_NAMES[cidx]
            result[f'param{i}'] = float(popt[i])
        return result
    except Exception:
        return None


def load_datasets():
    os.makedirs(INTERMEDIATE_DIR, exist_ok=True)
    df_pheno = pd.read_excel(PHYSIOLOGY_PATH, sheet_name='Sheet1')
    if 'sampleID' in df_pheno.columns:
        df_pheno = df_pheno.set_index('sampleID')
    df_pheno.index = df_pheno.index.astype(str)
    df_org = pd.read_excel(ORGANELLE_ALL_PATH, sheet_name='Sheet1')
    df_org = df_org.iloc[:, 1:]
    df_org = df_org.set_index(df_org.columns[0]).T
    df_org.index = df_org.index.astype(str)
    df_org.columns = df_org.columns.astype(str)
    train_mask = df_pheno['source'].astype(str).str.contains('Rosemary', case=False, na=False) & df_pheno.index.str.startswith('prot.')
    val_mask = df_pheno['source'].astype(str).str.contains('Jianye', case=False, na=False)
    df_pheno_train = df_pheno[train_mask]
    df_pheno_val = df_pheno[val_mask]
    df_org_train = df_org.loc[df_org.index.intersection(df_pheno_train.index)]
    df_org_val = df_org.loc[df_org.index.intersection(df_pheno_val.index)]
    common_cols = [c for c in df_org_train.columns if c in df_org_val.columns]
    df_org_train = df_org_train[common_cols]
    df_org_val = df_org_val[common_cols]
    return df_org_train, df_pheno_train, df_org_val, df_pheno_val, common_cols


def run_target(target, df_org_train, df_pheno_train, df_org_val, df_pheno_val, feature_names):
    if target == 'dilution rate (/h)':
        exclude_samples = [f'prot.{i}' for i in range(22, 43)]
        o_train = df_org_train[~df_org_train.index.isin(exclude_samples)]
        p_train = df_pheno_train[~df_pheno_train.index.isin(exclude_samples)]
        o_val = df_org_val[~df_org_val.index.isin(exclude_samples)]
        p_val = df_pheno_val[~df_pheno_val.index.isin(exclude_samples)]
    else:
        o_train = df_org_train
        p_train = df_pheno_train
        o_val = df_org_val
        p_val = df_pheno_val
    common_train = o_train.index.intersection(p_train.index)
    common_val = o_val.index.intersection(p_val.index)
    x_train_df = o_train.loc[common_train, feature_names]
    x_val_df = o_val.loc[common_val, feature_names]
    y_train = p_train.loc[common_train, target]
    y_val = p_val.loc[common_val, target]
    if target == 'qGlucose (mmol/gDW h)':
        y_val = y_val.abs()
    valid_cols = [c for c in feature_names if x_train_df[c].notna().all() and x_val_df[c].notna().all()]
    x_train_df = x_train_df[valid_cols]
    x_val_df = x_val_df[valid_cols]
    feature_names = list(x_train_df.columns)
    x_train = x_train_df.to_numpy(dtype=float)
    x_val = x_val_df.to_numpy(dtype=float)
    y_train = y_train.to_numpy(dtype=float)
    y_val = y_val.to_numpy(dtype=float)
    all_results = []
    for k in range(MIN_FEATURES, MAX_FEATURES + 1):
        total = math.comb(len(feature_names), k)
        heap = []
        combo_iter = combinations(range(len(feature_names)), k)
        with ProcessPoolExecutor(
            max_workers=MAX_WORKERS,
            initializer=_init_worker,
            initargs=(x_train, y_train, x_val, y_val, feature_names),
        ) as executor:
            for result in tqdm(executor.map(_fit_combo, combo_iter, chunksize=128), total=total, desc=f'{target} {k}特征', leave=False):
                if result is None:
                    continue
                score = result['Val_adj_R2']
                if np.isnan(score):
                    continue
                if len(heap) < TOP_K_PER_FEATURE_COUNT:
                    heapq.heappush(heap, (score, result))
                elif score > heap[0][0]:
                    heapq.heapreplace(heap, (score, result))
        top_results = [x[1] for x in sorted(heap, key=lambda t: t[0], reverse=True)]
        for item in top_results:
            item['目标变量'] = target
            item['特征数量'] = k
            all_results.append(item)
    return all_results


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
    df_org_train, df_pheno_train, df_org_val, df_pheno_val, feature_names = load_datasets()
    all_results = []
    for target in tqdm(TARGET_VARIABLES, desc='目标变量进度'):
        target_results = run_target(target, df_org_train, df_pheno_train, df_org_val, df_pheno_val, feature_names)
        all_results.extend(target_results)
    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values(['目标变量', 'Val_adj_R2'], ascending=[True, False]).reset_index(drop=True)
    for i in range(1, 9):
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
