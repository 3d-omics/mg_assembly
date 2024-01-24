#!/usr/bin/env python3

"""genomad_summary_to_gff.py: convert the _summary.tsv from geNomad to GFF3"""

import argparse

import pandas as pd


parser = argparse.ArgumentParser(
    prog="genomad_summary_to_gff",
    description="Convert geNomad's summary.tsv to GFF3 format",
)

parser.add_argument(
    "-i",
    "--input-summary",
    dest="summary",
    help="geNomad's *_summary.tsv",
    required=True,
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

    summary = pd.read_table(args.summary)

    summary["seqid"] = summary.seq_name
    summary["source"] = "genomad"
    summary["type"] = "CDS"
    summary["start"] = 1
    summary["end"] = summary.length
    if "virus_score" in summary.columns:
        summary["score"] = summary.virus_score
    elif "plastid_score" in summary.columns:
        summary["score"] = summary.plastid_score
    else:
        summary["score"] = 0
    summary["strand"] = "+"
    summary["phase"] = "."
    summary["attributes"] = [f"ID={x};Name={x}" for x in summary.seq_name]

    with open(args.gff, mode="w", encoding="utf-8") as gff:
        gff.write("##gff-version 3\n")

    summary[gff_colnames].to_csv(
        args.gff, sep="\t", header=False, mode="a", index=False
    )
