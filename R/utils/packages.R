# packages.R
# Package loading utilities

load_packages <- function(packages) {
  invisible(lapply(
    c(
      # core_packages
      "tidyverse",
      "data.table",
      "ggpubr",
      "ggfortify",
      "rstatix",
      "broom",
      "readxl",

      # model_packages
      "car",
      "nlme",
      "MASS",
      "MuMIn",
      "sjPlot",
      "class",
      "caret",

      # diversity_packages
      "vegan",
      "hillR",
      "geosphere",
      "indicspecies",

      # ggplot_extensions
      "ggtext",
      "ggpmisc",
      "MetBrewer",

      # table_packages
      "gt",
      "huxtable",
      "flextable",
      "broom.mixed",
      "officer",
      "knitr",

      # map_packages
      "ggmap",
      "ggspatial",

      # phyloseq_packages,
      "phyloseq",
      "microeco",
      "file2meco",
      "metagMisc",
      "microbiome",
      "mirlyn",

      # misc_packages
      "parallelly",
      "conflicted"
    ),
    function(pkg) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  ))
}


# Conflict preferences ----
set_conflicts <- function() {
  conflict_prefer("select", "dplyr")
  conflict_prefer("desc", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("alpha", "ggplot2")
}
