# setup.R
# Main orchestration script - sources modular utilities

# Source modular scripts ----
lapply(
  c(
    "packages.R",
    "load_data.R",
    "derived_data.R"
  ),
  function(utils) {
    source(here::here("R/utils", utils))
  }
)

# Load packages and set conflicts ----
load_packages()
set_conflicts()

# Custom functions ----
list.files(here::here("R/functions"), full.names = TRUE) |>
  purrr::walk(source)


# Comparison Lists ----
species_comparisons <- list(
  c("M. nasutus", "M. laciniatus"),
  c("M. nasutus", "M. guttatus"),
  c("M. laciniatus", "M. guttatus")
)

elevation_comparisons <- list(
  c("LOW", "MID"),
  c("LOW", "HIGH"),
  c("MID", "HIGH")
)
