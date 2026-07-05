## Data Format Requirements

All input data must follow the formats described below. Place your files in `data/` before running the pipeline.

---

### Genotype Data

PLINK2 format (`.pgen` / `.pvar` / `.psam`):

```
data/
├── your_data.pgen
├── your_data.pvar
└── your_data.psam
```

The `.psam` file must have an `#IID` column (no FID required):

```
#IID    SEX
HG00096  1
HG00097  2
```

---

### Phenotype File

Tab-separated, with a header row. Must contain `FID`, `IID`, and at least one phenotype column:

```
FID    IID        Expression
0      HG00096    -0.169414
0      HG00097    -0.152333
0      HG00099     0.909230
```

- `FID` can be all zeros
- `IID` must match the sample IDs in the `.psam` file
- Phenotype values should be continuous and ideally pre-normalized (mean=0, sd=1)

---

### Covariate File

Two supported formats:

**Long format** (rows = samples, preferred):

```
IID        PC1       PC2       sex
HG00096   -0.0240   0.0127    1
HG00097   -0.0243   0.0221    2
```

**Wide format** (rows = covariates, cols = samples): the pipeline will auto-detect and transpose this format.

```
           HG00096   HG00097   ...
PC1        -0.0240   -0.0243
PC2         0.0127    0.0221
sex         1         2
```

- `sex` should be encoded as `1` (male) / `2` (female)
- PC columns should be named `PC1`, `PC2`, etc.
- PEER/InferredCov columns are supported but not recommended for PRS analysis 

---

### HapMap3 Reference File

Required for LDpred2. Download `map_hm3_plus.rds` from the `bigsnpr` package resources and place it at:

```
data/map_hm3_plus.rds
```

This file provides the HapMap3 SNP list used for variant filtering. It must contain columns: `chr`, `pos`, `pos_hg38`, `a0`, `a1`, `rsid`.

---

### Notes

**Genome build**: all data must use **GRCh38 (hg38)**. The pipeline does not perform liftover.

**Sample IDs**: IID must be consistent across `.psam`, phenotype, and covariate files. Mismatches will cause samples to be silently dropped.

**Covariates**: the pipeline supports any combination of continuous and binary covariates. For standard GWAS and PRS analysis, we recommend using genetic PCs (`PC1`, `PC2`, ...) and `sex` only. Including expression-derived covariates (e.g. PEER/InferredCov) will reduce PRS predictive power.

**Multi-allelic variants and Indels**: the pipeline automatically filters out Indels via `--snps-only` in the GWAS step. Multi-allelic variants are rare in standard array/WGS data processed with standard QC pipelines.