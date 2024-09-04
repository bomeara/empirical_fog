library(targets)
library(tarchetypes)
library(crew)

source("functions.R")
source("_packages.R")


list(

  tar_target(geiger_data, get_data()),
  tar_target(all_data, get_Alencar_et_al_data(geiger_data)),
  tar_target(all_ouwie_results, run_datasets(all_data)),
  tar_target(all_geiger_results, run_datasets_geiger(all_data)),
  tar_target(all_ouwie_df, summarize_ouwie_results(all_ouwie_results)),
  tar_target(all_geiger_df, summarize_geiger_results(all_geiger_results)),
  tar_target(all_df, organize_final_df(dplyr::bind_rows(all_ouwie_df, all_geiger_df))),
  tar_target(all_df_csv, write.csv(all_df, "all_df.csv"))
)