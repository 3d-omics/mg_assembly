#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input-folder",
  type = "character",
  dest = "input_folder",
  help = "Folder containing the *.tsv files"
)

parser$add_argument(
  "-o", "--output-file",
  type = "character",
  dest = "output_file",
  help = "Output TSV file"
)

args <- parser$parse_args()
input_folder <- args$input_folder
output_file <- args$output_file
output_folder <- dirname(output_file)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

files <- list.files(args$input_folder, pattern = "*.tsv", full.names = TRUE)

# files <- list.files("results/metabin/coverm/contig/", pattern = "*.tsv", full.names = TRUE)

rename_cols <- function(x) {
  colnames(x)[1] <- "sequence_id"
  return(x)
}

nonempty_files <-
  files %>%
  map(function(x) read_tsv(x, col_types = cols()), .progress = TRUE) %>%
  keep(function(x) nrow(x) > 0)

if (length(nonempty_files) > 0) {
  nonempty_files %>%
    map(rename_cols) %>%
    map(function(x) pivot_longer(x, -sequence_id, names_to = "library", values_to = "counts")) %>%
    bind_rows() %>%
    mutate(library = str_split(library, " ") %>% map_chr(1)) %>%
    group_by(sequence_id, library) %>%
    summarise(counts = sum(counts)) %>%
    pivot_wider(names_from = "library", values_from = "counts", values_fill = NA) %>%
    write_tsv(output_file)
} else {
  write_tsv(x = tibble(Contig = NA), file = output_file)
}
