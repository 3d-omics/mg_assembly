include: "concoct_functions.smk"


rule concoct_cut_up_fasta_one:
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        assembly_10k=CONCOCT / "prepare" / "{assembly_id}.cut.fa",
        bed_10k=CONCOCT / "prepare" / "{assembly_id}.cut.bed",
    log:
        CONCOCT / "prepare" / "{assembly_id}.cut.log",
    conda:
        "../../envs/bin/concoct.yml"
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


rule concoct_coverage_table_one:
    input:
        bams=get_bams_for_concoct_binning,
        bais=get_bais_for_concoct_binning,
        bed_10k=CONCOCT / "prepare" / "{assembly_id}.cut.bed",
    output:
        coverage=CONCOCT / "prepare" / "{assembly_id}.coverage.tsv",
    log:
        CONCOCT / "prepare" / "{assembly_id}.coverage.log",
    conda:
        "../../envs/bin/concoct.yml"
    shell:
        """
        concoct_coverage_table.py \
            {input.bed_10k} \
            {input.bams} \
        > {output.coverage} \
        2>> {log}
        """


rule concoct_run_one:
    input:
        assembly_10k=CONCOCT / "prepare" / "{assembly_id}.cut.fa",
        coverage=CONCOCT / "prepare" / "{assembly_id}.coverage.tsv",
    output:
        out_dir=directory(CONCOCT / "run/{assembly_id}/"),
    log:
        CONCOCT / "run/{assembly_id}.log",
    conda:
        "../../envs/bin/concoct.yml"
    shell:
        """
        concoct \
            --composition_file {input.assembly_10k} \
            --coverage_file {input.coverage} \
            --basename {output.out_dir} \
        2>> {log} 1>&2
        """


rule concoct_merge_cutup_clustering_one:
    input:
        run_dir=CONCOCT / "run/{assembly_id}",
    output:
        clustering_merged=CONCOCT / "merge" / "{assembly_id}.csv",
    log:
        CONCOCT / "merge" / "{assembly_id}.log",
    conda:
        "../../envs/bin/concoct.yml"
    shell:
        """
        merge_cutup_clustering.py \
            {input.run_dir}/clustering_gt1000.csv \
        > {output.clustering_merged} \
        2>> {log} 1>&2
        """


rule concoct_extract_fasta_bins_one:
    input:
        assembly=ASSEMBLE_BOWTIE2 / "{assembly_id}.fa",
        clustering_merged=CONCOCT / "merge" / "{assembly_id}.csv",
    output:
        bins=CONCOCT / "fasta_bins" / "{assembly_id}/",
    log:
        CONCOCT / "fasta_bins/{assembly_id}.log",
    conda:
        "../../envs/bin/concoct.yml"
    shell:
        """
        extract_fasta_bins.py \
            {input.assembly} \
            {output.clustering_merged} \
            --output_path {output.bins} \
        2> {log} 1>&2
        """


rule concoct:
    input:
        [CONCOCT / "fasta_bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
