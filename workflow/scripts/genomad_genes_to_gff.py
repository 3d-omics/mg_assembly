#!/usr/bin/env python3

"""genomad_genes_to_gff.py: convert the _genes.tsv from geNomad to GFF3"""

import argparse

import pandas as pd


parser = argparse.ArgumentParser(
    prog="genomad_genes_to_gff",
    description="Convert geNomad's genes.tsv to GFF3 format",
)

parser.add_argument(
    "-i", "--input-genes", dest="genes", help="geNomad's *_genes.tsv", required=True
)

parser.add_argument(
    "-o",
    "--output-gff",
    dest="gff",
    help="GFF file with the gene coordinates",
    required=True,
)

if __name__ == "__main__":
    args = parser.parse_args()

    gff_colnames = "seqid source type start end score strand phase attributes".split(
        " "
    )

    genes = pd.read_table(args.genes)

    genes["seqid"] = genes.gene.str.rsplit("_", n=1).str[0]
    genes["source"] = "genomad"
    genes["type"] = "CDS"
    genes["score"] = genes.bitscore.fillna(".")
    genes["strand"] = ["+" if x == 1 else "-" for x in genes.strand]
    genes["phase"] = "."
    genes["attributes"] = [f"ID={x};Name={x}" for x in genes["gene"]]

    with open(args.gff, mode="w", encoding="utf-8") as gff:
        gff.write("##gff-version 3\n")

    genes[gff_colnames].to_csv(args.gff, sep="\t", header=False, mode="a", index=False)
