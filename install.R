# R dependencies installation script
packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "bigstatsr",
  "bigsnpr",
  "magrittr"
)

# check and install missing packages 
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

cat("Done!\n")