#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(bigsnpr)
library(ggplot2)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
(NCORES <- nb_cores())
library(magrittr)


# ── default parameters ────────────────────────────────────────
defaults <- list(
  data_dir     = "results/intermediate_data",
  gwas_file    = " ",
  hm3_dir      = "data/map_hm3_plus.rds",
  out_dir      = "data/results/LDPred2",
  val_prefix   = "val",
  test_prefix  = "test",
  covar_col    = " ",
  pheno_name   = " ",
  maf_thr      = 0.01,
  random_state = 42
)


# ── parse arguments ───────────────────────────────────────────
parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  params <- defaults  

  if (length(args) > 0) {
    if (length(args) %% 2 != 0) stop("Arguments must be in --key value pairs")
    keys <- gsub("^--", "", args[seq(1, length(args), 2)])
    vals <- args[seq(2, length(args), 2)]

    # check if parameter name is correct
    invalid <- setdiff(keys, names(defaults))
    if (length(invalid) > 0) {
      stop(paste("Unknown argument(s):", paste(invalid, collapse = ", "),
                 "\nValid arguments:", paste(names(defaults), collapse = ", ")))
    }

    for (i in seq_along(keys)) {
      params[[keys[i]]] <- vals[i]
    }
  }

  params$maf_thr      <- as.numeric(params$maf_thr)
  params$random_state <- as.integer(params$random_state)

  return(params)
}
params <- parse_args(defaults)


# ── validate required files ───────────────────────────────────
check_file <- function(path, name) {
  if (!file.exists(path)) stop(paste(name, "not found:", path))
}

data_dir     <- params$data_dir
gwas_dir     <- file.path(data_dir, params$gwas_file)
hap_map3_dir <- params$hm3_dir
out_dir      <- params$out_dir
val_prefix   <- file.path(data_dir, params$val_prefix)
test_prefix  <- file.path(data_dir, params$test_prefix)

check_file(gwas_dir,                       "--gwas_file")
check_file(hap_map3_dir,                   "--hm3_dir")
check_file(paste0(val_prefix,  ".bed"),    "--val_prefix")
check_file(paste0(test_prefix, ".bed"),    "--test_prefix")
check_file(paste0(val_prefix,  "_pheno.txt"), "val pheno")
check_file(paste0(test_prefix, "_pheno.txt"), "test pheno")

val_pheno_dir  <- paste0(val_prefix,  "_pheno.txt")
val_covar_dir  <- paste0(val_prefix,  "_covar.txt")
test_pheno_dir <- paste0(test_prefix, "_pheno.txt")
test_covar_dir <- paste0(test_prefix, "_covar.txt")
covar_col      <- strsplit(params$covar_col, ",")[[1]]
pheno_name     <- params$pheno_name
cols_need      <- c("IID", pheno_name, covar_col)
maf_thr        <- params$maf_thr
random_state   <- params$random_state

cat("Parameters:\n")
for (k in names(params)) cat(" ", k, "=", params[[k]], "\n")



# ── load data ─────────────────────────────────────────────────
## GWAS summary statistics
gwas <- fread(gwas_dir)
rsID_available <- any(grepl("^rs", gwas$ID))

## Hap Map3
info <- readRDS(hap_map3_dir)

## reference panel 
if (!file.exists(paste0(data_dir, "/base.bed"))) {
  cmd <- paste0("./plink2 --pfile ", data_dir, "/base --make-bed --out ", data_dir, "/base")
  system(cmd)
}
if (!file.exists(paste0(data_dir, "/base.rds"))){
  snp_readBed(paste0(data_dir, "/base.bed"))
}
obj.bigSNP <- snp_attach(paste0(data_dir, "/base.rds"))
genotype <- obj.bigSNP$genotypes
chromosome <- unique(obj.bigSNP$map$chr)


## val set
if (!file.exists(paste0(val_prefix, ".rds"))) {
  snp_readBed(paste0(val_prefix, ".bed"))
}
val_set <- snp_attach(paste0(val_prefix, ".rds"))
val_G <- val_set$genotypes
val_pheno <- read.table(val_pheno_dir, header = TRUE)
val_covar <- read.table(val_covar_dir, header = TRUE)
val_data <- left_join(val_pheno, val_covar, by = "IID")[, cols_need]

## test set
if (!file.exists(paste0(test_prefix, ".rds"))) {
  snp_readBed(paste0(test_prefix, ".bed"))
}
test_set <- snp_attach(paste0(test_prefix, ".rds"))
test_G <- test_set$genotypes
test_pheno <- read.table(test_pheno_dir, header = TRUE)
test_covar <- read.table(test_covar_dir, header = TRUE)
test_data <- left_join(test_pheno, test_covar, by = "IID")[, cols_need]

# ── preprocess ────────────────────────────────────────────────
## filter out Hap Map3 variants 
sumstats <- data.table(
  chr     = gwas$`#CHROM`,
  pos     = gwas$POS,
  rsid    = gwas$ID,
  a1      = gwas$A1,        # effect allele
  a0      = gwas$REF,       # reference allele
  n_eff   = gwas$OBS_CT[1], # GWAS sample size
  beta_se = gwas$SE,
  p       = gwas$P,
  beta    = gwas$BETA
)
if (rsID_available) {
  sumstats <- sumstats[sumstats$rsid %in% info$rsid, ]
  cat("rsID available in GWAS data!\n")
  cat("SNPs in sumstats after HapMap3 filter:", nrow(sumstats), "\n")
} else {
  sumstats <- sumstats[sumstats$pos %in% info$pos_hg38, ]
  cat("rsID not found in GWAS data, use POS instead.\n")
  cat("SNPs in sumstats after HapMap3 filter:", nrow(sumstats), "\n")
}

## match referance panel & summary stats, use SNPs in summary stats only
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
if (rsID_available) {
  df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)
} else {
  df_beta <- snp_match(sumstats, map)
}

## save what we shrank for Ridge Model
write.csv(df_beta$rsid, paste0(data_dir, "/snps_after_shrink.csv"), row.names=FALSE)


# ── calculate LD matrix ───────────────────────────────────────
## filter out variants with small MAF
ind.row <- rows_along(genotype)
maf <- snp_MAF(genotype, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
df_beta <- df_beta[maf > maf_thr, ]
cat("After MAF QC,", nrow(df_beta), " variants remain.\n")
## create and open a temporary file
tmp_dir <- paste0(data_dir, "/tmp-data")
if (dir.exists(tmp_dir)) {
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}
dir.create(tmp_dir)
tmp <- tempfile(tmpdir = tmp_dir)
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

## convert physical positions (in bp) to genetic positions (in cM)
POS2 <- snp_asGeneticPos(map$chr, map$pos, dir = tmp_dir, ncores = NCORES)

## LD calculation
cat("Start calculating LD matrix")
for (i in seq_along(chromosome)) {
  chr <- chromosome[i]
  cat("Processing chr", chr, "\n")
  ind.chr <- which(df_beta$chr == chr)     ## indices in 'df_beta'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]  ## indices in 'genotype'
  
  corr0 <- snp_cor(genotype, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (i == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
cat("Done! Size of LD matrix:", file.size(corr$sbk) / 1024^3,
    "GB, with dimension = ", paste(dim(corr), collapse = " x "), "\n")


# ── grid search ───────────────────────────────────────────────
## get estimated h2
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
if (h2_est < 0){
  h2_est <- 0.1 # h2_est should be positive, if sample generates a negative one, we replace it with a small value
}

## make grid
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
grid_param <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

## train the effect weights
set.seed(random_state)  # reproducible
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid_param, ncores = NCORES)


# ── evaluation ───────────────────────────────────────────────
## make PRS estimation
val_pred_grid <- big_prodMat(val_G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
y_val <- val_data[[pheno_name]]

y_test <- test_data[[pheno_name]]
test_pred_grid <- big_prodMat(test_G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

## make baseline prediction
base_mod <- lm(as.formula(paste(pheno_name, "~ .")), data = val_data[, -which(names(val_data) == "IID")])
val_base_pred  <- predict(base_mod)
val_base_rmse  <- sqrt(mean((y_val - val_base_pred)^2))
test_base_pred <- predict(base_mod, newdata = test_data[, -which(names(test_data) %in% c("IID", pheno_name))])
test_base_rmse <- sqrt(mean((y_test - test_base_pred)^2))


## evaluation loop
val_corr_list  <- c()
val_rmse_list  <- c()
test_corr_list <- c()
test_rmse_list <- c()



for (i in 1:ncol(beta_grid)) {
  if (all(is.na(val_pred_grid[, i]))) {
    val_corr_list  <- c(val_corr_list,  NA)
    val_rmse_list  <- c(val_rmse_list,  NA)
    test_corr_list <- c(test_corr_list, NA)
    test_rmse_list <- c(test_rmse_list, NA)
    next
  }
  
  # validation: fit calibration model
  val_df <- val_data[, -which(names(val_data) == "IID")]
  val_df$PRS <- val_pred_grid[, i]
  mod <- lm(as.formula(paste(pheno_name, "~ .")), data = val_df)
  
  val_pred   <- predict(mod)
  val_corr   <- cor(y_val, val_df$PRS)
  val_rmse   <- sqrt(mean((y_val - val_pred)^2))
  val_corr_list <- c(val_corr_list, val_corr)
  val_rmse_list <- c(val_rmse_list, val_rmse)
  
  # test: predict using calibration model
  test_df <- test_data[, -which(names(test_data) %in% c("IID", pheno_name))]
  test_df$PRS <- test_pred_grid[, i]
  test_pred  <- predict(mod, newdata = test_df)
  test_corr  <- cor(y_test, test_df$PRS)
  test_rmse  <- sqrt(mean((y_test - test_pred)^2))
  test_corr_list <- c(test_corr_list, test_corr)
  test_rmse_list <- c(test_rmse_list, test_rmse)
}

best_idx <- which.max(abs(val_corr_list))

# ── plot ──────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_plot <- function(plt_data, y_col, y_lab, title, best_idx, sub_title, baseline = NULL) {
  p <- ggplot(plt_data, aes(x = p, y = .data[[y_col]], color = as.factor(h2))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    geom_point(data = plt_data[best_idx, ],
               color = "red", size = 4, shape = 4) +
    scale_x_log10(breaks = 10^(-5:0), minor_breaks = plt_data$p) +
    facet_wrap(~ sparse, labeller = label_both) +
    labs(title = title,
         subtitle = paste0(sub_title, " = ",
                           round(plt_data[[y_col]][best_idx], 4)),
         y = y_lab, x = "Causal Variant Fraction (p)", color = "h2") +
    theme(legend.position = "top", panel.spacing = unit(1, "lines"))
  
  if (!is.null(baseline)) {
    p <- p + geom_hline(yintercept = baseline, linetype = "dashed",
                        color = "black", linewidth = 0.8) +
      annotate("text", x = min(plt_data$p), y = baseline,
               label = paste0("baseline=", round(baseline, 4)),
               vjust = -0.5, hjust = 0, size = 3)
  }
  return(p)
}

val_plt_data  <- cbind(grid_param, corr = val_corr_list,  rmse = val_rmse_list)
test_plt_data <- cbind(grid_param, corr = test_corr_list, rmse = test_rmse_list)

plots <- list(
  val_corr  = make_plot(val_plt_data,  "corr", "Correlation(Pheno & PRS)", "Validation - Correlation", 
                        best_idx, "best correlation"),
  
  val_rmse  = make_plot(val_plt_data,  "rmse", "RMSE",                     "Validation - RMSE",        
                        best_idx, "best RMSE from validation", baseline = val_base_rmse),
  
  test_corr = make_plot(test_plt_data, "corr", "Correlation(Pheno & PRS)", "Test - Correlation",       
                        best_idx, "best correlation from validation"),
  
  test_rmse = make_plot(test_plt_data, "rmse", "RMSE",                     "Test - RMSE",              
                        best_idx, "best RMSE from validation", baseline = test_base_rmse)
)

for (name in names(plots)) {
  ggsave(paste0(out_dir, "/", name, ".png"), plot = plots[[name]],
         width = 30, height = 15, units = "cm", dpi = 300)
}