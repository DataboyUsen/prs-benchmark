"""Data splitting script for PRS benchmarking pipeline."""

import pandas as pd
import numpy as np
import argparse
import subprocess
import os

def split_data(prefix, data_dir, out_dir, pheno_file, covar_file,
               base_ratio, val_ratio, test_ratio, random_state):

    assert abs(base_ratio + val_ratio + test_ratio - 1.0) < 1e-6, \
        "Ratios must sum to 1"
    os.makedirs(out_dir, exist_ok=True)

    # Read sample list
    psam = pd.read_csv(os.path.join(data_dir, f"{prefix}.psam"), sep="\t")
    psam_shuffled = psam.sample(frac=1, random_state=random_state).reset_index(drop=True)
    n = len(psam_shuffled)

    n_base  = int(n * base_ratio)
    n_val = int(n * val_ratio)

    splits = {
        "base":  psam_shuffled.iloc[:n_base],
        "val": psam_shuffled.iloc[n_base:n_base+n_val],
        "test":  psam_shuffled.iloc[n_base+n_val:]
    }
    print(f"Total: {n} | Base: {len(splits['base'])} | "
          f"Validation: {len(splits['val'])} | Test: {len(splits['test'])}")

    # 保存样本列表并切割 pgen
    for name, subset in splits.items():
        iid_col = "#IID" if "#IID" in psam.columns else psam.columns[0]
        sample_file = os.path.join(out_dir, f"{name}_samples.txt")
        subset[[iid_col]].to_csv(sample_file, index=False, header=False)

        pfile_out = os.path.join(out_dir, name)
        pfile_in = os.path.join(data_dir, prefix)
        if name == "base":
            cmd = f"./plink2 --pfile {pfile_in} --keep {sample_file} --make-pgen --out {pfile_out}"  # base data use PLINK2 style for GWAS computing 
        else:
            cmd = f"./plink2 --pfile {pfile_in} --keep {sample_file} --make-bed --out {pfile_out}"   # other data use PLINK1 style for compatibility
        subprocess.run(cmd, shell=True, check=True)
        print(f"{name} pgen created: {pfile_out}")

    # 切割 pheno
    if pheno_file:
        pheno = pd.read_csv(os.path.join(data_dir, pheno_file), sep="\t")
        iid_col = "IID" if "IID" in pheno.columns else pheno.columns[0]
        for name, subset in splits.items():
            ids = subset["#IID"].values if "#IID" in psam.columns else subset.iloc[:, 0].values
            out = pheno[pheno[iid_col].isin(ids)]
            out.to_csv(os.path.join(out_dir, f"{name}_pheno.txt"), sep="\t", index=False)
        print("Pheno files created")

    # 切割 covar
    if covar_file:
        covar = pd.read_csv(os.path.join(data_dir, covar_file), sep="\t")
        # 判断是否需要转置（宽格式）
        iid_col = "IID" if "IID" in covar.columns else None
        if iid_col is None:
            # 宽格式，样本在列，转置
            covar = covar.set_index(covar.columns[0]).T
            covar.index.name = "IID"
            covar = covar.reset_index()
            iid_col = "IID"
        for name, subset in splits.items():
            ids = subset["#IID"].values if "#IID" in psam.columns else subset.iloc[:, 0].values
            out = covar[covar[iid_col].isin(ids)]
            out.to_csv(os.path.join(out_dir, f"{name}_covar.txt"), sep="\t", index=False)
        print("Covar files created")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split PLINK2 data into base/val/test")
    parser.add_argument("--prefix",       type=str,   required=True,  help="PLINK2 file prefix")
    parser.add_argument("--data-dir",     type=str,   default=".",    help="Input data directory")
    parser.add_argument("--out-dir",      type=str,   default=".",    help="Output directory")
    parser.add_argument("--pheno",        type=str,   default=None,   help="Phenotype file")
    parser.add_argument("--covar",        type=str,   default=None,   help="Covariate file")
    parser.add_argument("--base-ratio",   type=float, default=0.6)
    parser.add_argument("--val-ratio",    type=float, default=0.2)
    parser.add_argument("--test-ratio",   type=float, default=0.2)
    parser.add_argument("--random-state", type=int,   default=42)
    args = parser.parse_args()

    split_data(
        prefix=args.prefix,
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        pheno_file=args.pheno,
        covar_file=args.covar,
        base_ratio=args.base_ratio,
        val_ratio=args.val_ratio,
        test_ratio=args.test_ratio,
        random_state=args.random_state
    )