#!/usr/bin/env python

"""
Merge contig-to-bin TSV files from different binning tools (CONCOCT, MaxBin2, MetaBAT2).
"""

import argparse
import gzip
import os
import sys

# pylint: disable=import-error
from Bio import SeqIO


def get_bin_id(file_path):
    """Get the bin ID from the file path"""
    basename = os.path.basename(file_path)
    if basename.endswith(".gz"):
        basename = basename[:-3]
    for ext in [".fa", ".fasta", ".fna"]:
        if basename.endswith(ext):
            return basename[: -len(ext)]
    return basename


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Merge contig to bin files from different binners."
    )
    parser.add_argument("--concoct", help="CONCOCT bin directory", required=True)
    parser.add_argument("--maxbin2", help="MaxBin2 bin directory", required=True)
    parser.add_argument("--metabat2", help="MetaBAT2 bin directory", required=True)
    parser.add_argument("--output", help="Output TSV file", required=True)
    return parser.parse_args()


def validate_directories(binners):
    """Validate that all provided input paths are directories"""
    for binner_dir in binners.values():
        if not os.path.isdir(binner_dir):
            sys.exit(f"Error: {binner_dir} is not a directory or does not exist.")


def extract_contigs_from_fasta(file_path):
    """Generator yielding contig IDs from a FASTA file (can handle .gz)"""
    opener = gzip.open if file_path.endswith(".gz") else open
    with opener(file_path, "rt", encoding="utf-8") as f:
        for record in SeqIO.parse(f, "fasta"):
            yield record.id


def process_binner_directory(binner_dir, binner_name, out_f):
    """List files in a binner directory and write contig-to-bin mappings to out_f"""

    fasta_extensions = (
        ".fa",
        ".fasta",
        ".fna",
        ".fna",
        ".fa.gz",
        ".fasta.gz",
        ".fna.gz",
        ".faa.gz",
    )

    for file_name in sorted(os.listdir(binner_dir)):
        file_path = os.path.join(binner_dir, file_name)
        if os.path.isfile(file_path) and file_name.endswith(fasta_extensions):
            bin_id = get_bin_id(file_path)
            for contig_id in extract_contigs_from_fasta(file_path):
                out_f.write(f"bin_{bin_id}\t{contig_id}\t{binner_name}\n")


def main():
    """Main execution function"""
    args = parse_args()

    binners = {
        "concoct": args.concoct,
        "maxbin2": args.maxbin2,
        "metabat2": args.metabat2,
    }

    validate_directories(binners)

    with open(args.output, "w", encoding="utf-8") as out_f:
        for binner_name, binner_dir in binners.items():
            process_binner_directory(binner_dir, binner_name, out_f)


if __name__ == "__main__":
    main()
