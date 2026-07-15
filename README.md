# PRS-Benchmark

Benchmarking pipeline for Polygenic Risk Score (PRS) methods: C+T (PRSice-2), LDpred2, and Ridge Regression.

## Requirements

**Software**
- Python ≥ 3.9
- R ≥ 4.0
- PLINK2
- PRSice-2 (Mac binary + PRSice.R)

**Python packages**
```bash
pip install -r requirements.txt
```

**R packages**
```r
source("install.R")
```

---

## Data

Place your data files in `data/`:
```
data/
├── your_data.pgen
├── your_data.pvar
├── your_data.psam
├── pheno.txt        # columns: IID, <pheno_name>
└── covar.txt        # wide format: rows=covariates, cols=samples
```

---

## Pipeline

### Step 0: Split Data

Split PLINK2 genotype, phenotype, and covariate files into base / train / test sets.
- **base**: PLINK2 format, used for GWAS
- **val**: PLINK1 format, used for hyper-parameter selection and PRS calibration
- **test**: PLINK1 format, used for final evaluation

```bash
python pipeline/0_split_data.py \
  --prefix <your_data> \
  --data-dir data/ \
  --out-dir results/intermediate_data \
  --pheno <pheno_file>.txt \
  --covar <covar_file>.txt \
  --base-ratio 0.6 \
  --val-ratio 0.2 \
  --test-ratio 0.2
```

| Parameter | Description |
|---|---|
| `--prefix` | PLINK2 file prefix |
| `--data-dir` | Directory containing input data |
| `--out-dir` | Directory for output files |
| `--pheno` | Phenotype filename, optional |
| `--covar` | Covariate filename, optional. Supports both wide format (rows=covariates, cols=samples) and long format (rows=samples) |
| `--base-ratio` | Proportion for GWAS base set (default: 0.6) |
| `--val-ratio` | Proportion for hyper-parameter selection and PRS calibration (default: 0.2) |
| `--test-ratio` | Proportion for final evaluation (default: 0.2) |
| `--random-state` | Random seed (default: 42) |

Output:
```
results/intermediate_data/base.pgen / .pvar / .psam
results/intermediate_data/val.bed / .bim / .fam
results/intermediate_data/test.bed  / .bim / .fam
results/intermediate_data/base_pheno.txt / val_pheno.txt / test_pheno.txt
results/intermediate_data/base_covar.txt / val_covar.txt / test_covar.txt
```

---

### Step 1: GWAS

Run genome-wide association study on the base set using PLINK2.

```bash
bash pipeline/1_gwas.sh \
  --pfile results/intermediate_data/base \
  --pheno results/intermediate_data/base_pheno.txt \
  --pheno-name <pheno_name> \
  --covar results/intermediate_data/base_covar.txt \
  --covar-cols <covariate_list> \
  --out results/intermediate_data/gwas_result
```

| Parameter | Description |
|---|---|
| `--pfile` | PLINK2 base genotype prefix |
| `--pheno` | Phenotype file |
| `--pheno-name` | Column name of the phenotype to analyse |
| `--covar` | Covariate file, optional |
| `--covar-cols` | Comma-separated covariate column names to include (e.g. `PC1,PC2,sex`) |
| `--out` | Output prefix |
| `--threads` | Number of threads (default: 16) |

Output: `results/intermediate_data/gwas_result.<pheno_name>.glm.linear`

---

### Step 2: C+T (PRSice-2)

Place `PRSice_mac` and `PRSice.R` in `pipeline/` before running.

```bash
bash pipeline/2_ct_prsice.sh \
  --prsice PRSice_mac \
  --base-gwas results/intermediate_data/gwas_result.<pheno_name>.glm.linear \
  --target-val results/intermediate_data/val \
  --target-test results/intermediate_data/test \
  --pheno-val results/intermediate_data/val_pheno.txt \
  --pheno-test results/intermediate_data/test_pheno.txt \
  --pheno-col <pheno_name> \
  --bar-levels "5e-8,1e-6,1e-4,1e-3,0.01,0.05,0.1,0.2,0.5,1" \
  --out results/CT/ct_result
```

| Parameter | Description |
|---|---|
| `--prsice` | Path to PRSice-2 binary (e.g. `PRSice_mac`) |
| `--base-gwas` | GWAS summary statistics file |
| `--target-val` | PLINK1 prefix for validation set |
| `--target-test` | PLINK1 prefix for test set |
| `--pheno-val` | Phenotype file for validation set |
| `--pheno-test` | Phenotype file for test set |
| `--pheno-col` | Column name of the phenotype |
| `--bar-levels` | Comma-separated p-value thresholds to test |
| `--base-maf` | Minimum MAF filter for base GWAS SNPs (default: 0.01) |
| `--out` | Output prefix |

Output: `results/CT/ct_result_val.all_score`, `results//CT/ct_result_test.all_score`

---

### Step 2b: C+T Evaluation

Evaluate C+T PRS scores generated in Step 2. Fits a calibration model on the validation set and evaluates RMSE and correlation on both validation and test sets across all p-value thresholds.

```bash
Rscript pipeline/2b_ct_eval.R \
  --data_dir results/intermediate_data \
  --out_dir results/CT \
  --val_prefix ct_result_val \
  --test_prefix ct_result_test \
  --pheno_name <pheno_name> \
  --covar_col <covariate_list>
```

| Parameter | Description | Default |
|---|---|---|
| `--data_dir` | Directory containing pheno and covar files
| `--out_dir` | Directory containing CT PRS score files and for saving plots 
| `--val_prefix` | Filename prefix of validation set score file (default: `ct_result_val`)
| `--test_prefix` | Filename prefix of test set score file (default: `ct_result_test`)
| `--pheno_name` | Column name of the phenotype 
| `--covar_col` | Comma-separated covariate column names (e.g. `PC1,PC2,sex`)

Output:
```
results/CT/val_corr.png    # correlation across thresholds, validation set
results/CT/val_rmse.png    # RMSE across thresholds, validation set
results/CT/test_corr.png   # correlation across thresholds, test set
results/CT/test_rmse.png   # RMSE across thresholds, test set
```

---

### Step 3: LDpred2 and Lassosum2

Runs LDpred2-grid on the base set genotype and GWAS summary statistics. Automatically detects whether rsID is available in the GWAS data; if not, falls back to position-based matching with HapMap3. Evaluates PRS performance on both validation and test sets. At the same time, a lassosum2-grid (defult 120 parameter combinations) will also be trained on base set. 

```bash
Rscript pipeline/3_ldpred2_lassosum.R \
  --data_dir results/intermediate_data \
  --gwas_file gwas_result.<pheno_name>.glm.linear \
  --hm3_dir data/map_hm3_plus.rds \
  --out_dir results \
  --val_prefix val \
  --test_prefix test \
  --pheno_name <pheno_name> \
  --covar_col <covariate_list> \
  --maf_thr 0.01 \
  --random_state 42
```

| Parameter | Description | Default |
|---|---|---|
| `--data_dir` | Directory containing base/val/test plink files and pheno/covar files | `results/intermediate_data` |
| `--gwas_file` | GWAS summary statistics filename (inside `--data_dir`) | ` ` |
| `--hm3_dir` | Path to HapMap3 RDS file (`map_hm3_plus.rds`) | `data/map_hm3_plus.rds` |
| `--out_dir` | Directory for output plots | `results/LDPred2` |
| `--val_prefix` | Filename prefix for validation set (`.bed`, `_pheno.txt`, `_covar.txt`) | `val` |
| `--test_prefix` | Filename prefix for test set (`.bed`, `_pheno.txt`, `_covar.txt`) | `test` |
| `--pheno_name` | Column name of the phenotype | ` ` |
| `--covar_col` | Comma-separated covariate column names (e.g. `PC1,PC2,sex`) | ` ` |
| `--maf_thr` | MAF threshold for filtering variants | `0.01` |
| `--random_state` | Random seed for reproducibility | `42` |

**Note:** The script expects the following files in `--data_dir`:
```
base.pgen / .pvar / .psam    # GWAS reference panel (PLINK2)
base.bed  / .bim  / .fam    # auto-generated if not present
val.bed   / .bim  / .fam
val_pheno.txt / val_covar.txt
test.bed  / .bim  / .fam
test_pheno.txt / test_covar.txt
```

**Note:** `map_hm3_plus.rds` is required for HapMap3 SNP filtering. It can be downloaded from the `bigsnpr` package resources.

Output:
```
results/LDPred2/val_corr.png     # correlation across grid, validation set
results/LDPred2/val_rmse.png     # RMSE across grid, validation set
results/LDPred2/test_corr.png    # correlation across grid, test set
results/LDPred2/test_rmse.png    # RMSE across grid, test set
results/lassosum/val_corr.png    # correlation across grid, validation set
results/lassosum/val_rmse.png    # RMSE across grid, validation set
results/lassosum/test_corr.png   # correlation across grid, test set
results/lassosum/test_rmse.png   # RMSE across grid, test set
results/intermediate_data/snps_after_shrink.csv      # SNP list after HapMap3
```

---

### Step 4: Ridge

Runs 2 Ridge models with full variants and shrunk variants as previous method `LDPred-2` used. Manually set whether to use validation set; if `--include_val=True`, combines validation set and base set as training data; if not, only use base set **(currently only support `--include_val=True`)**. Evaluates PRS performance on test sets.

```bash
python pipeline/4_ridge.py \
  --data_dir results/intermediate_data \
  --out_dir results/Ridge \
  --use_snp results/intermediate_data/snps_after_shrink.csv \
  --covar_cols <covariate_list> \
  --pheno_name <pheno_name> \
  --include_val True \
  --alphas "0.1,1,10,100,1000"
```

| Parameter | Description | Default |
|---|---|---|
| `--data_dir` | Directory containing base/val/test plink files and pheno/covar files | `results/intermediate_data` |
| `--out_dir` | Directory for output plots | `results/Ridge` |
| `--use_snp` | Path to SNPs used after HapMap3-shrinkage, the program would select variants in this file | `results/intermediate_data/snps_after_shrink.csv` |
| `--covar_cols` | Comma-separated covariate column names (e.g. `PC1,PC2,sex`) | ` ` |
| `--pheno_name` | Column name of the phenotype | ` ` |
| `--include_val` | If include validation set as training data | `True` |
| `--alphas` | Comma-separated alphas, (e.g. `0.1,1,10`)  | If empty, uses `np.logspace(-1,10,30)`

**Note:** `include_val=False` is currently not supported.

Output:
```
results/Ridge/ridge_test_rmse_all_variants.png     # RMSE across grid when all variants are used, test set
results/Ridge/ridge_test_rmse_hapmap3.png          # RMSE across grid when HapMap3 shrunk variants are used, test set
```

---




