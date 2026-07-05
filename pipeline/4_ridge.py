#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import os

# ── parse arguments ───────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(description="Ridge PRS benchmarking")
    parser.add_argument("--data_dir",    type=str, default="results/intermediate_data")
    parser.add_argument("--out_dir",     type=str, default="results/Ridge")
    parser.add_argument("--use_snp",     type=str, default="results/intermediate_data/snps_after_shrink.csv")
    parser.add_argument("--covar_cols",  type=str, default="")
    parser.add_argument("--pheno_name",  type=str, default="")
    parser.add_argument("--include_val", type=str, default="True",  help="Include val set in training (True/False)")
    parser.add_argument("--alphas",      type=str, default="",       help="Comma-separated alphas, e.g. '0.1,1,10'. If empty, uses np.logspace(-1,10,30)")
    return parser.parse_args()

args        = parse_args()
data_dir    = args.data_dir
out_dir     = args.out_dir
use_snp     = args.use_snp
covar_cols  = args.covar_cols.split(",")
pheno_name  = args.pheno_name
include_val = args.include_val.lower() == "true"
alphas      = np.array([float(a) for a in args.alphas.split(",")]) if args.alphas else np.logspace(-1, 10, 30)

os.makedirs(out_dir, exist_ok=True)
print(f"include_val: {include_val}")
print(f"alphas: {alphas}")

# ── helper functions ──────────────────────────────────────────
def load_plink_data(prefix, pheno_file, covar_file, covar_cols, pheno_name):
    """Load plink1 bed + pheno + covar, align by IID, return merged DataFrame."""
    G = read_plink1_bin(f"{prefix}.bed", verbose=False)
    G_df = pd.DataFrame(
        G.values,
        columns=G.snp.values
    )
    G_df.insert(0, "IID", [s.decode() if isinstance(s, bytes) else s for s in G.iid.values])

    pheno = pd.read_csv(pheno_file, sep="\t")
    covar = pd.read_csv(covar_file, sep="\t")[["IID"] + covar_cols]

    # merge on IID to ensure alignment
    df = G_df.merge(pheno[["IID", pheno_name]], on="IID") \
             .merge(covar, on="IID")

    # reorder: IID | genotype | covariates | phenotype
    geno_cols = list(G.snp.values)
    df = df[["IID"] + geno_cols + covar_cols + [pheno_name]]
    return df

def extract_Xy(df, covar_cols, pheno_name):
    """Extract X (genotype + covariate) and y from DataFrame."""
    geno_cols = [c for c in df.columns if c not in ["IID"] + covar_cols + [pheno_name]]
    X = df[geno_cols + covar_cols].values.astype(np.float32)
    y = df[pheno_name].values.astype(np.float32)
    # fill missing genotype with column mean
    col_means = np.nanmean(X, axis=0)
    for j in range(X.shape[1]):
        X[np.isnan(X[:, j]), j] = col_means[j]
    return X, y

def run_ridge(X_train, y_train, X_test, y_test, alphas, out_dir, tag=""):
    """Train Ridge for each alpha, evaluate on test, plot RMSE."""
    rmse_list = []
    for alpha in alphas:
        model = Ridge(alpha=alpha)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        rmse_list.append(rmse)
        print(f"  alpha={alpha:.2e}  RMSE={rmse:.4f}")

    best_idx  = np.argmin(rmse_list)
    best_rmse = rmse_list[best_idx]
    print(f"Best alpha: {alphas[best_idx]:.2e}, RMSE: {best_rmse:.4f}")

    # baseline: covariate only
    n_geno   = X_train.shape[1] - len(covar_cols)
    X_cov_tr = X_train[:, n_geno:]
    X_cov_te = X_test[:,  n_geno:]
    from sklearn.linear_model import LinearRegression
    lr = LinearRegression().fit(X_cov_tr, y_train)
    baseline_rmse = np.sqrt(mean_squared_error(y_test, lr.predict(X_cov_te)))
    print(f"Baseline RMSE: {baseline_rmse:.4f}")

    # plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(np.log10(alphas), rmse_list, color="steelblue", linewidth=2, marker="o", markersize=3)
    ax.scatter(np.log10(alphas[best_idx]), best_rmse, color="red", s=80, zorder=5,
               label=f"best alpha={alphas[best_idx]:.2e}, RMSE={best_rmse:.4f}")
    ax.axhline(baseline_rmse, color="black", linestyle="--",
               label=f"baseline={baseline_rmse:.4f}")
    ax.set_xlabel("log10(alpha)")
    ax.set_ylabel("RMSE")
    ax.set_title(f"Ridge PRS - Test RMSE {tag}")
    ax.legend()
    plt.tight_layout()
    fname = os.path.join(out_dir, f"ridge_test_rmse{tag}.png")
    plt.savefig(fname, dpi=300)
    plt.close()
    print(f"Plot saved: {fname}")

# ── round 1: all variants ─────────────────────────────────────
print("\n=== Round 1: All variants ===")
df_base = load_plink_data(
    os.path.join(data_dir, "base"),
    os.path.join(data_dir, "base_pheno.txt"),
    os.path.join(data_dir, "base_covar.txt"),
    covar_cols, pheno_name
)
df_val = load_plink_data(
    os.path.join(data_dir, "val"),
    os.path.join(data_dir, "val_pheno.txt"),
    os.path.join(data_dir, "val_covar.txt"),
    covar_cols, pheno_name
)
df_test = load_plink_data(
    os.path.join(data_dir, "test"),
    os.path.join(data_dir, "test_pheno.txt"),
    os.path.join(data_dir, "test_covar.txt"),
    covar_cols, pheno_name
)

print(f"df_base:  {df_base.shape}")
print(f"df_val:   {df_val.shape}")
print(f"df_test:  {df_test.shape}")

# assemble train
df_train = pd.concat([df_base, df_val], ignore_index=True) if include_val else df_base

X_train, y_train = extract_Xy(df_train, covar_cols, pheno_name)
X_test,  y_test  = extract_Xy(df_test,  covar_cols, pheno_name)

if include_val:
    run_ridge(X_train, y_train, X_test, y_test, alphas, out_dir, tag="_all_variants")

# ── round 2: HapMap3 variants ─────────────────────────────────
print("\n=== Round 2: HapMap3 variants ===")
snp_list = pd.read_csv(use_snp)["x"].values
print(f"HapMap3 SNPs to use: {len(snp_list)}")

def filter_snps(df, snp_list, covar_cols, pheno_name):
    """Keep only HapMap3 SNP columns."""
    geno_cols = [c for c in df.columns if c not in ["IID"] + covar_cols + [pheno_name]]
    keep_geno = [c for c in geno_cols if c in snp_list]
    print(f"  Matched {len(keep_geno)} / {len(snp_list)} HapMap3 SNPs")
    return df[["IID"] + keep_geno + covar_cols + [pheno_name]]

df_base = filter_snps(df_base, snp_list, covar_cols, pheno_name)
df_val  = filter_snps(df_val,  snp_list, covar_cols, pheno_name)
df_test = filter_snps(df_test, snp_list, covar_cols, pheno_name)

print(f"df_base after shrink:  {df_base.shape}")
print(f"df_val  after shrink:  {df_val.shape}")
print(f"df_test after shrink:  {df_test.shape}")

df_train = pd.concat([df_base, df_val], ignore_index=True) if include_val else df_base

X_train, y_train = extract_Xy(df_train, covar_cols, pheno_name)
X_test,  y_test  = extract_Xy(df_test,  covar_cols, pheno_name)

if include_val:
    run_ridge(X_train, y_train, X_test, y_test, alphas, out_dir, tag="_hapmap3")