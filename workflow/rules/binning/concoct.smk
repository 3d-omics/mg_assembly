include: "concoct_functions.smk"


rule binning_concoct_cut_up_fasta_one:
    input:
        assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
    output:
        assembly_10k=CONCOCT / "cut" / "{assembly_id}.fa",
        bed_10k=CONCOCT / "cut" / "{assembly_id}.bed",
    log:
        CONCOCT / "cut" / "{assembly_id}.log",
    conda:
        "../../envs/binning/concoct.yml"
    shell:
        """
        cut_up_fasta.py \
            {input.assembly} \
            --chunk_size 10000 \
            --overlap_size 0 \
            --merge_last \
            --bedfile {output.bed_10k} \
        > {output.assembly_10k} \
        2> {log}
        """


rule binning_concoct_coverage_table_one:
    input:
        bams=get_bams_for_concoct_binning,
        bais=get_bais_for_concoct_binning,
        bed_10k=CONCOCT / "cut" / "{assembly_id}.bed",
    output:
        coverage=CONCOCT / "coverage" / "{assembly_id}/coverage_table.tsv",
    log:
        CONCOCT / "coverage" / "{assembly_id}.log",
    conda:
        "../../envs/binning/concoct.yml"
    shell:
        """
        concoct_coverage_table.py \
            {input.bed_10k} \
            {input.bams} \
        > {output.coverage} \
        2>> {log}
        """


rule binning_concoct_run_one:
    input:
        assembly_10k=CONCOCT / "cut" / "{assembly_id}.fa",
        coverage=CONCOCT / "coverage" / "{assembly_id}/coverage_table.tsv",
    output:
        out_dir=directory(CONCOCT / "{assembly_id}/concoct_output"),
    log:
        CONCOCT / "{assembly_id}/concoct.log",
    conda:
        "../../envs/binning/concoct.yml"
    shell:
        """
        echo $(which concoct) 2>> {log} 1>&2
        concoct \
            --composition_file {input.assembly_10k} \
            --coverage_file {input.coverage} \
            --basename {output.out_dir} \
        2>> {log} 1>&2
        """


rule binning_concoct_merge_cutup_clustering_one:
    input:
        out_dir=CONCOCT / "{assembly_id}/concoct_output",
    output:
        clustering_merged=CONCOCT / "merge" / "{assembly_id}/clustering_merged.csv",
    log:
        CONCOCT / "merge" / "{assembly_id}.log",
    conda:
        "../../envs/binning/concoct.yml"
    shell:
        """
        merge_cutup_clustering.py \
            {output.outdir}/clustering_gt1000.csv \
        > {output.clustering_merged} \
        2>> {log} 1>&2
        """


rule binning_concoct_extract_fasta_bins_one:
    input:
        assembly=MEGAHIT_RENAMING / "{assembly_id}.fa",
        clustering_merged=CONCOCT / "merge" / "{assembly_id}/clustering_merged.csv",
    output:
        bins=CONCOCT / "fasta_bins" / "{assembly_id}/",
    log:
        CONCOCT / "fasta_bins/{assembly_id}.log",
    conda:
        "../../envs/binning/concoct.yml"
    shell:
        """
        extract_fasta_bins.py \
            {input.assembly} \
            {output.clustering_merged} \
            --output_path {output.bins} \
        2> {log} 1>&2
        """


rule binning_concoct_all:
    input:
        [CONCOCT / "fasta_bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
