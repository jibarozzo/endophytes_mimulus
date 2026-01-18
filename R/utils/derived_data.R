#' Prepare Derived Data
#'
#' Creates species-specific sample name vectors and other derived objects
#' from the loaded datasets.
#'
#' @param data Named list from `load_datasets()`.
#' @return Named list with original data plus derived objects.
#' @export
prepare_derived_data <- function(data) {
  data$leaf_traits |>
    filter(!Species == "M. bicolor")
  data$leaf_traits_noB <- data$leaf_traits |>
    filter(!Species == "M. bicolor")

  data$final_names_methodB <- colnames(data$asv_avgdist)
  data$asv_matrix <- otu_table(data$ps_clean_3) |> # ASV matrix
    as.data.frame() |>
    select(contains(data$final_names_methodB)) |> # Samples from method B randomization
    as.matrix()

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
