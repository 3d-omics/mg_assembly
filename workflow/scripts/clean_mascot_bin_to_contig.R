library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input-file",
  type = "character",
  dest = "input_file",
  help = "Bin to contig file from magscot"
)

parser$add_argument(
  "-o", "--output-file",
  type = "character",
  dest = "output_file",
  help = "clean bin to contig file"
)

args <- parser$parse_args()
input_file <- args$input_file
output_file <- args$output_file
output_folder <- dirname(output_file)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

raw_magscot <- read_tsv(input_file)

raw_magscot %>%
  mutate(
    binnew = binnew %>%
      str_split("/") %>%
      map_chr(-1) %>%
      str_remove("magscot_cleanbin_"),
  ) %>%
  separate(
    col = contig,
    into = c("assembly_id", "contig_id"),
    remove = FALSE
  ) %>%
  mutate(
    seqname = str_glue("{assembly_id}.bin_{binnew}.{contig_id}")
  ) %>%
  write_tsv(output_file)
