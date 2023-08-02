#!/usr/bin/env bash

set -euo pipefail

output_folder="resources"

wget \
    --continue \
    https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz


tar \
    -xvzf \
    gtdbtk_r214_data.tar.gz \
    -C $output_folder
