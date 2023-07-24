
library(Nonpareil)
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Input and/or output files are missing", call. = FALSE)
}

Nonpareil.set(args[1], plot = FALSE) %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(., var = "sample") %>%
  write.table(
    args[2],
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
