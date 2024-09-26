#!/usr/bin/Rscript

library(tidyverse)
library(argparse)

parser <- argparse::ArgumentParser()

parser$add_argument(
  "-i", "--input-folder",
  type = "character",
  dest = "input_folder",
  help = "Folder containing the *.npo files"
)

parser$add_argument(
  "-t", "--output-tsv",
  type = "character",
  dest = "output_tsv",
  help = "Output TSV file"
)



args <- parser$parse_args()
input_folder <- args$input_folder
output_tsv <- args$output_tsv

dir.create(output_tsv %>% dirname(), showWarnings = FALSE, recursive = TRUE)

nonempty_files <-
  list.files(args$input_folder, pattern = "*.npo", full.names = TRUE) %>%
  data.frame(file = .) %>%
  mutate(size = file.info(file)$size) %>%
  filter(size > 0) %>%
  pull(file)

if (length(nonempty_files) > 0) {

  nonempty_files %>%
    Nonpareil::Nonpareil.set(plot = FALSE) %>%
    summary() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    write_tsv(output_tsv)

} else {

  write_tsv(data.frame(), output_tsv)

}
