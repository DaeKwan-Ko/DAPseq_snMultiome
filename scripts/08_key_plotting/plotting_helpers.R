# plotting_helpers.R
# Helper functions shared by key plotting scripts.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

theme_manuscript <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

save_plot <- function(plot, filename, width = 4, height = 3) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
  message("Saved plot: ", filename)
}

read_optional_tsv <- function(path) {
  if (file.exists(path)) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    message("Missing optional input: ", path)
    tibble()
  }
}

standardize_clade <- function(x) {
  dplyr::case_when(
    x %in% c("A", "ECA") ~ "ECA",
    x %in% c("B", "ECB") ~ "ECB",
    TRUE ~ as.character(x)
  )
}

standardize_condition <- function(x) {
  dplyr::case_when(
    x %in% c("Control", "Ctrl", "Ctl", "control") ~ "Control",
    x %in% c("Salinity", "NaCl", "Salt", "salt", "salinity") ~ "Salinity",
    TRUE ~ as.character(x)
  )
}
