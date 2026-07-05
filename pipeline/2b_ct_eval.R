#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(bigstatsr)

# ── default parameters ────────────────────────────────────────
defaults <- list(
  data_dir    = "results/intermediate_data",
  out_dir     = "results/CT",
  val_prefix  = "ct_result_val",
  test_prefix = "ct_result_test",
  pheno_name  = "Expression",
  covar_col   = "PC1,PC2,PC3,PC4,PC5,sex"
)

# ── parse arguments ───────────────────────────────────────────
parse_args <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  params <- defaults
  if (length(args) > 0) {
    if (length(args) %% 2 != 0) stop("Arguments must be in --key value pairs")
    keys <- gsub("^--", "", args[seq(1, length(args), 2)])
    vals <- args[seq(2, length(args), 2)]
    invalid <- setdiff(keys, names(defaults))
    if (length(invalid) > 0) {
      stop(paste("Unknown argument(s):", paste(invalid, collapse = ", "),
                 "\nValid arguments:", paste(names(defaults), collapse = ", ")))
    }
    for (i in seq_along(keys)) params[[keys[i]]] <- vals[i]
  }
  return(params)
}

params      <- parse_args(defaults)
data_dir    <- params$data_dir
out_dir     <- params$out_dir
pheno_name  <- params$pheno_name
covar_col   <- strsplit(params$covar_col, ",")[[1]]
cols_need   <- c("IID", pheno_name, covar_col)

# ── load data ─────────────────────────────────────────────────
val_scores  <- fread(file.path(out_dir, paste0(params$val_prefix,  ".all_score")))
test_scores <- fread(file.path(out_dir, paste0(params$test_prefix, ".all_score")))

val_pheno  <- read.table(file.path(data_dir, "val_pheno.txt"),  header = TRUE)
val_covar  <- read.table(file.path(data_dir, "val_covar.txt"),  header = TRUE)
test_pheno <- read.table(file.path(data_dir, "test_pheno.txt"),   header = TRUE)
test_covar <- read.table(file.path(data_dir, "test_covar.txt"),   header = TRUE)

val_data  <- left_join(val_pheno,  val_covar,  by = "IID")[, cols_need]
test_data <- left_join(test_pheno, test_covar, by = "IID")[, cols_need]

# align sample order
val_data  <- val_data[match(val_scores$IID,  val_data$IID), ]
test_data <- test_data[match(test_scores$IID, test_data$IID), ]

y_val  <- val_data[[pheno_name]]
y_test <- test_data[[pheno_name]]

# ── get threshold columns ─────────────────────────────────────
thresh_cols <- grep("^Pt_", names(val_scores), value = TRUE)

# ── baseline ──────────────────────────────────────────────────
val_base_df <- val_data[, -which(names(val_data) == "IID")]
base_mod    <- lm(as.formula(paste(pheno_name, "~ .")), data = val_base_df)

val_base_rmse  <- sqrt(mean((y_val  - predict(base_mod))^2))
test_base_df   <- test_data[, -which(names(test_data) %in% c("IID", pheno_name))]
test_base_rmse <- sqrt(mean((y_test - predict(base_mod, newdata = test_base_df))^2))

cat("Baseline RMSE (val):",  round(val_base_rmse,  4), "\n")
cat("Baseline RMSE (test):", round(test_base_rmse, 4), "\n")

# ── evaluation loop ───────────────────────────────────────────
val_corr_list  <- c()
val_rmse_list  <- c()
test_corr_list <- c()
test_rmse_list <- c()

for (thresh in thresh_cols) {
  if (all(is.na(val_scores[[thresh]]))) {
    val_corr_list  <- c(val_corr_list,  NA)
    val_rmse_list  <- c(val_rmse_list,  NA)
    test_corr_list <- c(test_corr_list, NA)
    test_rmse_list <- c(test_rmse_list, NA)
    next
  }
  
  # validation: fit calibration model
  val_df      <- val_data[, -which(names(val_data) == "IID")]
  val_df$PRS  <- val_scores[[thresh]]
  mod         <- lm(as.formula(paste(pheno_name, "~ .")), data = val_df)
  val_pred    <- predict(mod)
  val_corr_list <- c(val_corr_list, cor(y_val, val_df$PRS))
  val_rmse_list <- c(val_rmse_list, sqrt(mean((y_val - val_pred)^2)))
  
  # test: predict using calibration model
  test_df      <- test_data[, -which(names(test_data) %in% c("IID", pheno_name))]
  test_df$PRS  <- test_scores[[thresh]]
  test_pred    <- predict(mod, newdata = test_df)
  test_corr_list <- c(test_corr_list, cor(y_test, test_df$PRS))
  test_rmse_list <- c(test_rmse_list, sqrt(mean((y_test - test_pred)^2)))
}

best_idx <- which.max(abs(val_corr_list))
cat("Best threshold:", thresh_cols[best_idx], "\n")
cat("Val  corr:", round(val_corr_list[best_idx],  4), "\n")
cat("Val  RMSE:", round(val_rmse_list[best_idx],  4), "\n")
cat("Test corr:", round(test_corr_list[best_idx], 4), "\n")
cat("Test RMSE:", round(test_rmse_list[best_idx], 4), "\n")

# ── plot ──────────────────────────────────────────────────────
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# p values from threshold column names
p_vals <- as.numeric(gsub("Pt_", "", thresh_cols))

make_plot <- function(corr_list, rmse_list, p_vals, best_idx,
                      base_rmse, set_name) {
  plt_data <- data.frame(p = p_vals, corr = corr_list, rmse = rmse_list)
  
  p_corr <- ggplot(plt_data, aes(x = p, y = corr)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_point(data = plt_data[best_idx, ], color = "red", size = 4, shape = 4) +
    scale_x_log10() +
    labs(title = paste0("C+T - ", set_name, " Correlation"),
         subtitle = paste0("Validation gives corr = ", round(corr_list[best_idx], 4)),
         x = "P-value Threshold", y = "Correlation (PRS vs Phenotype)")
  
  p_rmse <- ggplot(plt_data, aes(x = p, y = rmse)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_point(data = plt_data[best_idx, ], color = "red", size = 4, shape = 4) +
    geom_hline(yintercept = base_rmse, linetype = "dashed", color = "black") +
    annotate("text", x = min(p_vals), y = base_rmse,
             label = paste0("baseline=", round(base_rmse, 4)),
             vjust = -0.5, hjust = 0, size = 3) +
    scale_x_log10() +
    labs(title = paste0("C+T - ", set_name, " RMSE"),
         subtitle = paste0("Validation gives RMSE = ",  round(rmse_list[best_idx], 4)),
         x = "P-value Threshold", y = "RMSE")
  
  list(corr = p_corr, rmse = p_rmse)
}

val_plots  <- make_plot(val_corr_list,  val_rmse_list,  p_vals, best_idx, val_base_rmse,  "Validation")
test_plots <- make_plot(test_corr_list, test_rmse_list, p_vals, best_idx, test_base_rmse, "Test")

ggsave(file.path(out_dir, "val_corr.png"),  val_plots$corr,  width = 20, height = 10, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "val_rmse.png"),  val_plots$rmse,  width = 20, height = 10, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "test_corr.png"), test_plots$corr, width = 20, height = 10, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "test_rmse.png"), test_plots$rmse, width = 20, height = 10, units = "cm", dpi = 300)
cat("Plots saved to", out_dir, "\n")
