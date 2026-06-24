from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import Ridge
from sklearn.feature_selection import RFE
from sklearn.metrics import r2_score, mean_squared_error
import numpy as np
import pandas as pd
import os
from datetime import datetime
import multiprocessing as mp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, '..', '..', 'data')
ORIGINAL_DATA_DIR = os.path.join(DATA_DIR, 'original_data')
INTERMEDIATE_DIR = os.path.join(DATA_DIR, 'intermediate')
organelle_path = os.path.join(ORIGINAL_DATA_DIR, "main_organelle_fraction_260414.xlsx")
genes_path = os.path.join(ORIGINAL_DATA_DIR, "df_genes.csv")
physiology_path = os.path.join(ORIGINAL_DATA_DIR, "physiology_collection.xlsx")

organelles_main = ["mitochondrion", "nucleus", "endoplasmic reticulum",
 "fungal-type vacuole",  "plasma membrane", "ribosome", "cytosol"]

X_full = None
y_full = None
out_dir = None
header = "fold\tmin_protein_count\tR2_cv\tRMSE_cv\tR2_test\tRMSE_test\tproteins\n"

def run_seed(split_seed):
    output_txt_file = os.path.join(out_dir, f"10fold_R2_results_seed{split_seed}.txt")
    X_train_cv, X_test, y_train_cv, y_test = train_test_split(X_full, y_full, test_size=0.1, random_state=split_seed)

    kf = KFold(n_splits=10, shuffle=True, random_state=split_seed)
    splits = list(kf.split(X_train_cv))

    with open(output_txt_file, "w", encoding="utf-8") as f:
        f.write(header)

    for fold_idx, (train_idx, val_idx) in enumerate(splits, 1):
        X_train, X_val = X_train_cv.iloc[train_idx], X_train_cv.iloc[val_idx]
        y_train, y_val = y_train_cv.iloc[train_idx], y_train_cv.iloc[val_idx]
        for n in range(1, 61):
            model = Ridge(random_state=split_seed)
            selector = RFE(model, n_features_to_select=n)
            selector.fit(X_train, y_train)
            selected_cols = X_train.columns[selector.support_]
            X_train_sel = X_train[selected_cols]
            X_val_sel = X_val[selected_cols]
            X_test_sel = X_test[selected_cols]
            model.fit(X_train_sel, y_train)
            y_pred_val = model.predict(X_val_sel)
            y_pred_test = model.predict(X_test_sel)
            r2_cv = r2_score(y_val, y_pred_val)
            rmse_cv = np.sqrt(mean_squared_error(y_val, y_pred_val))
            r2_test = r2_score(y_test, y_pred_test)
            rmse_test = np.sqrt(mean_squared_error(y_test, y_pred_test))
            proteins = list(selected_cols)
            txt_line = f"{fold_idx}\t{n}\t{r2_cv}\t{rmse_cv}\t{r2_test}\t{rmse_test}\t{';'.join(proteins)}\n"
            with open(output_txt_file, "a", encoding="utf-8") as f:
                f.write(txt_line)

def main():
    global X_full, y_full, out_dir
    seeds = [1, 2, 3, 4, 5]
    out_dir = os.path.join(INTERMEDIATE_DIR, "multi_seed")
    df_organelle_raw = pd.read_excel(organelle_path)
    df_organelle = df_organelle_raw[df_organelle_raw["compartment"].isin(organelles_main)].copy()
    df_organelle = df_organelle.drop(columns=[c for c in ["Unnamed: 0"] if c in df_organelle.columns])
    df_organelle = df_organelle.set_index("compartment").loc[organelles_main].T
    sample_names = df_organelle.index.astype(str).str.strip()
    df_genes = pd.read_csv(genes_path)
    df_physiology = pd.read_excel(physiology_path)
    df_sheet5 = pd.read_excel(physiology_path, sheet_name="Sheet5")
    if "condition_unique" not in df_sheet5.columns:
        raise ValueError("Sheet5 must contain condition_unique column")
    sheet5_conditions = (
        df_sheet5["condition_unique"]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .tolist()
    )

    if len(df_genes) != len(sample_names):
        raise ValueError(f"Sample count mismatch: df_genes={len(df_genes)}, df_organelle={len(sample_names)}")

    required_cols = {"sampleID", "condition_unique"}
    if not required_cols.issubset(df_physiology.columns):
        raise ValueError("260417_physiology_collection.xlsx must contain sampleID and condition_unique columns")

    df_physiology["sampleID"] = df_physiology["sampleID"].astype(str).str.strip()
    df_physiology["condition_unique"] = df_physiology["condition_unique"].astype(str).str.strip()
    condition_map = (
        df_physiology[["sampleID", "condition_unique"]]
        .dropna(subset=["sampleID", "condition_unique"])
        .drop_duplicates(subset=["sampleID"], keep="first")
        .set_index("sampleID")["condition_unique"]
    )
    condition_index = pd.Index([condition_map.get(s, np.nan) for s in sample_names], name="condition_unique")
    if pd.isna(condition_index).any():
        missing_ids = [sample_names[i] for i, v in enumerate(condition_index) if pd.isna(v)]
        raise ValueError(f"Missing condition_unique mapping for sampleID values: {missing_ids[:20]}")

    X_full = df_genes.copy()
    X_full.index = sample_names
    X_full = X_full.groupby(condition_index, sort=False).mean()
    X_full = X_full.reindex(sheet5_conditions)

    y_full = df_organelle.copy()
    y_full.index = sample_names
    y_full = y_full.groupby(condition_index, sort=False).mean()
    y_full = y_full.reindex(sheet5_conditions)

    invalid_mask = X_full.isna().all(axis=1) | y_full.isna().all(axis=1)
    if bool(invalid_mask.any()):
        invalid_conditions = X_full.index[invalid_mask].astype(str).tolist()
        raise ValueError(f"Some Sheet5 conditions have no data in X/y after grouping: {invalid_conditions[:30]}")

    os.makedirs(out_dir, exist_ok=True)

    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass

    procs = []
    for split_seed in seeds:
        p = mp.Process(target=run_seed, args=(split_seed,))
        p.start()
        procs.append(p)

    for p in procs:
        p.join()

    print("out_dir:", out_dir)


if __name__ == "__main__":
    main()
