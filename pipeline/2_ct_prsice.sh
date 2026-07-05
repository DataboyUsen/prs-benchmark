# ── default parameters ────────────────────────────────────────
PRSICE=""
BASE_GWAS=""
TARGET_VAL=""       # validation set plink prefix
TARGET_TEST=""      # test set plink prefix
PHENO_VAL=""
PHENO_TEST=""
PHENO_COL=""
BAR_LEVELS="5e-8,1e-6,1e-4,1e-3,0.01,0.05,0.1,0.2,0.5,1"
BASE_MAF="0.01"
OUT="ct_result"

# ── parse arguments ───────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --prsice)         PRSICE="$2";         shift 2 ;;
        --base-gwas)      BASE_GWAS="$2";      shift 2 ;;
        --target-val)     TARGET_VAL="$2";     shift 2 ;;
        --target-test)    TARGET_TEST="$2";    shift 2 ;;
        --pheno-val)      PHENO_VAL="$2";      shift 2 ;;
        --pheno-test)     PHENO_TEST="$2";     shift 2 ;;
        --pheno-col)      PHENO_COL="$2";      shift 2 ;;
        --bar-levels)     BAR_LEVELS="$2";     shift 2 ;;
        --base-maf)       BASE_MAF="$2";       shift 2 ;;
        --out)            OUT="$2";            shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── validate ──────────────────────────────────────────────────
if [[ -z "$PRSICE" || -z "$BASE_GWAS" || -z "$TARGET_VAL" || -z "$TARGET_TEST" ]]; then
    echo "Error: --prsice, --base-gwas, --target-val, --target-test are required"
    exit 1
fi

# ── run ───────────────────────────────────────────────────────
mkdir -p $(dirname $OUT)

run_prsice() {
    local TARGET=$1
    local PHENO=$2
    local OUT_PREFIX=$3

    Rscript PRSice.R \
        --prsice $PRSICE \
        --base $BASE_GWAS \
        --target $TARGET \
        --binary-target F \
        --pheno $PHENO \
        --pheno-col $PHENO_COL \
        --base-maf A1_FREQ:$BASE_MAF \
        --stat BETA \
        --beta \
        --snp ID \
        --A1 A1 \
        --pvalue P \
        --bar-levels $BAR_LEVELS \
        --fastscore \
        --no-regress \
        --all-score \
        --ignore-fid \
        --out $OUT_PREFIX
}

echo "Running C+T on validation set..."
run_prsice $TARGET_VAL $PHENO_VAL ${OUT}_val

echo "Running C+T on test set..."
run_prsice $TARGET_TEST $PHENO_TEST ${OUT}_test

echo "Done: ${OUT}_val.all_score and ${OUT}_test.all_score"