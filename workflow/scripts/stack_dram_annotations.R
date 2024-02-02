library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-o", "--output-file",
  type = "character",
  dest = "output_file",
  help = "output tsv file"
)

parser$add_argument(
  'input_files',
  metavar='INPUT',
  type="character",
  nargs='+',
  help='TSV files to join'
)

args <- parser$parse_args()




input_files <- args$input_files
output_file <- args$output_file

map(
  .x = input_files,
  .f = function(x) {
    read_tsv(
      file =x,
      col_types="cccccccccccccccccccccccccccccccccccccccc")
  }
) %>%
bind_rows() %>%
write_tsv(output_file)

