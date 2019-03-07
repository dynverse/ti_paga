#!/usr/local/bin/Rscript

# generate dataset with certain seed
set.seed(1)
data <- dyntoy::generate_dataset(
  id = "specific_example/paga",
  num_cells = 99,
  num_features = 101,
  model = "tree"
)
params <- list()

# add method specific args (if needed)
data$params <- list()

# write example dataset to file
file <- commandArgs(trailingOnly = TRUE)[[1]]
dyncli::write_h5(data[c("counts", "expression", "params", "prior_information")], file)
