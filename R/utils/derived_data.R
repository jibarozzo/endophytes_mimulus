#' Prepare Derived Data
#'
#' Creates species-specific sample name vectors and other derived objects
#' from the loaded datasets.
#'
#' @param data Named list from `load_datasets()`.
#' @return Named list with original data plus derived objects.
#' @export
prepare_derived_data <- function(data) {
  # Species-specific sample names (first 157 observations have ASV data)
  data$lac_names <- data$plant_traits |>
    dplyr::slice(1:157) |>
    dplyr::filter(Species == "M. laciniatus") |>
    dplyr::pull(Unique_ID)

  data$gut_names <- data$plant_traits |>
    dplyr::slice(1:157) |>
    dplyr::filter(Species == "M. guttatus") |>
    dplyr::pull(Unique_ID)

  data$nas_names <- data$plant_traits |>
    dplyr::slice(1:157) |>
    dplyr::filter(Species == "M. nasutus") |>
    dplyr::pull(Unique_ID)

  # ASV sample names from phyloseq object
  data$names_list <- colnames(data$ps_clean_3@otu_table)

  data
}
