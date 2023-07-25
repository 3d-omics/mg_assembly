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

nonempty_files <-
  files %>%
  map(function(x) read_tsv(x, col_types = cols()), .progress = TRUE) %>%
  keep(function(x) nrow(x) > 0)

if (length(nonempty_files) > 0) {
  nonempty_files %>%
  reduce(left_join) %>%
  write_tsv(output_file)
} else {
  write_tsv(x = tibble(Contig = NA), file = output_file)
}
