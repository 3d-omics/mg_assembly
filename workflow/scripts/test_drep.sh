#!/usr/bin/env bash
set -euo pipefail

mkdir results/prokaryotes/cluster/drep/separated_bins --parents
cp \
    resources/mock_assembly.fa.gz \
    results/prokaryotes/cluster/drep/separated_bins/bin1.fa.gz

cp \
    resources/mock_assembly.fa.gz \
    results/prokaryotes/cluster/drep/separated_bins/bin2.fa.gz

./run -j 8 --rerun-triggers mtime -- assemble__drep
