#!/usr/bin/env Rscript

parser <- argparse::ArgumentParser()

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
print(args)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# input_file <- "results/prokaryotes/cluster/magscot/all/magscot.refined.contig_to_bin.out"

raw_magscot <- readr::read_tsv(input_file)

bin_location <-
  raw_magscot$binnew[1] |>
  stringr::str_split("/", simplify = TRUE) |>
  purrr::map_chr(-1) |>
  length()

raw_magscot |>
  dplyr::mutate(
    binnew = binnew |>
      stringr::str_split("/") |>
      purrr::map_chr(bin_location) |>
      stringr::str_remove("magscot_cleanbin_")
  ) |>
  tidyr::separate(
    col = contig,
    into = c("tmp", "contig_id"),
    sep = "@",
    remove = FALSE
  ) |>
  tidyr::separate(
    col = tmp,
    into = c("assembly_id", "bin_id"),
    sep = ":",
    remove = TRUE
  ) |>
  dplyr::mutate(
    seqname = stringr::str_glue("{assembly_id}:bin_{binnew}@{contig_id}")
  ) |>
  readr::write_tsv(output_file)
