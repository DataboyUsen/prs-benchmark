#!/bin/bash
# pipeline/1_gwas.sh
# Usage: bash 1_gwas.sh --pfile base --pheno base_pheno.txt --covar base_covar.txt --out gwas_result

# ── default parameters ────────────────────────────────────────
PFILE=""
PHENO=""
PHENO_NAME=""        # e.g. "BMI", "Height"
COVAR=""
COVAR_COLS=""        # e.g. "PC1,PC2,sex"
OUT="gwas_result"
THREADS=16

# ── parse arguments ───────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --pfile)       PFILE="$2";       shift 2 ;;
        --pheno)       PHENO="$2";       shift 2 ;;
        --pheno-name)  PHENO_NAME="$2";  shift 2 ;;
        --covar)       COVAR="$2";       shift 2 ;;
        --covar-cols)  COVAR_COLS="$2";  shift 2 ;;
        --out)         OUT="$2";         shift 2 ;;
        --threads)     THREADS="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── validate ──────────────────────────────────────────────────
if [[ -z "$PFILE" || -z "$PHENO" ]]; then
    echo "Error: --pfile and --pheno are required"
    exit 1
fi

# ── build command ─────────────────────────────────────────────
CMD="./plink2 \
  --pfile $PFILE \
  --snps-only \
  --pheno $PHENO \
  --pheno-name $PHENO_NAME \
  --glm hide-covar \
  --threads $THREADS \
  --out $OUT"

# ── covariate----- ────────────────────────────────────────────
if [[ -n "$COVAR" ]]; then
    CMD="$CMD --covar $COVAR --covar-variance-standardize"
fi

if [[ -n "$COVAR_COLS" ]]; then
    CMD="$CMD --covar-name $COVAR_COLS"
fi

# ── run ───────────────────────────────────────────────────────
echo "Running GWAS..."
echo "$CMD"
eval $CMD
echo "Done: ${OUT}.${PHENO_NAME}.glm.linear"