#' Load Master Datasets
#'
#' Loads all cleaned and processed data files for statistical analyses.
#'
#' # Main data sets
#' The data cleaning and wrangling of `leaf_traits` and `plant_traits` data sets are in the `mim2_statistics.qmd` notebook, chunk `cleaning_chaping_data`. The data sets are saved as RDS files and CSV files in the `statistics` folder. The data cleaning and wrangling for phyloseq objects are in the `mim2_bioinformatics.qmd` notebook. The phyloseq objects are saved as RDS files and CSV files in the `taxonomy` folder.

#' The FEF community data is only available for 160/320 samples. For the downstream analyses, we eliminated *M. bicolor* or species "B" in the data set. The removal of this is due to a small sample size for leaf trait measurements and even smaller for FEF community data (*n* = 3). The removal is not explicit in `plant_traits` it is only explicit in the subset `leaf_traits_noB` and latter on in the `asv_matrix` that gives way for community analyses in [Community Diversity] in the mim2_statistics.qmd. We include a custom function for p-value formatting.
#'
#' @param path Character. Base path to the data directory.
#' @return Named list containing all loaded datasets.
#' @export
load_datasets <- function(path) {
  # RDS files to load
  rds_files <- list(
    leaf_traits = "data/clean_data/statistics/leaf_traits.rds",
    plant_traits = "data/clean_data/statistics/plant_traits.rds",
    plant_traits_MB = "data/clean_data/statistics/plant_traits_MB.rds",
    ps_clean_3 = "data/clean_data/taxonomy/02-TAXA_8450_phyloseq_nonsingletons_noB.rds",
    pseq_rrfb = "data/clean_data/statistics/pseq_rrfb.rds",
    rarefied_phyloseq = "data/clean_data/statistics/rarefied_phyloseq.rds",
    rrfy_hell_matrix = "data/clean_data/statistics/rrfy_hell_matrix.rds",
    asv_avgdist = "data/clean_data/statistics/asv_avgdist.rds",
    geo_distm = "data/clean_data/statistics/geo_distm.rds"
  )

  # CSV files to load
  csv_files <- list(
    ps_clean_3_df = "data/clean_data/taxonomy/02-TAXA_8450_phyloseq_nonsingletons_noB.csv"
  )

  # Load RDS files using purrr
  datasets <- purrr::map(rds_files, ~ readRDS(file.path(path, .x)))

  # Load CSV files
  csv_data <- purrr::map(csv_files, ~ read.csv(file.path(path, .x)))

  # Combine all datasets
  c(datasets, csv_data)
}
