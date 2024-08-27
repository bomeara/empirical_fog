library(targets)
library(tarchetypes)
library(crew)

source("functions.R")
source("_packages.R")


list(

  tar_target(geiger_data, get_data()),
  tar_target(all_results, run_datasets(geiger_data)),
  tar_target(all_df, summarize_results(all_results)),
  tar_target(all_df_csv, write.csv(all_df, "all_df.csv"))
)