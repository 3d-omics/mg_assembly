rule _assemble__concoct__cut_up_fasta:
    """Run concoct cut_up_fasta.py on one assembly."""
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
    output:
        assembly_10k=CONCOCT / "prepare" / "{assembly_id}.cut.fa",
        bed_10k=CONCOCT / "prepare" / "{assembly_id}.cut.bed",
    log:
        CONCOCT / "prepare" / "{assembly_id}.cut.log",
    conda:
        "concoct.yml"
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


rule _assemble__concoct__coverage_table:
    """Run concoct concoct_coverage_table.py on one assembly."""
    input:
        bams=get_bams_from_assembly_id,
        bais=get_bais_from_assembly_id,
        bed_10k=CONCOCT / "prepare" / "{assembly_id}.cut.bed",
    output:
        coverage=CONCOCT / "prepare" / "{assembly_id}.coverage.tsv",
    log:
        CONCOCT / "prepare" / "{assembly_id}.coverage.log",
    conda:
        "concoct.yml"
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        concoct_coverage_table.py \
            {input.bed_10k} \
            {input.bams} \
        > {output.coverage} \
        2>> {log}
        """


rule _assemble__concoct__run:
    """Run concoct on one assembly."""
    input:
        assembly_10k=CONCOCT / "prepare" / "{assembly_id}.cut.fa",
        coverage=CONCOCT / "prepare" / "{assembly_id}.coverage.tsv",
    output:
        args=CONCOCT / "run" / "{assembly_id}_args.txt",
        clustering=CONCOCT / "run" / "{assembly_id}_clustering_gt1000.csv",
        log=CONCOCT / "run" / "{assembly_id}_log.txt",
        original_data=CONCOCT / "run" / "{assembly_id}_original_data_gt1000.csv",
        components=CONCOCT / "run" / "{assembly_id}_PCA_components_data_gt1000.csv",
        transformed_data=CONCOCT
        / "run"
        / "{assembly_id}_PCA_transformed_data_gt1000.csv",
    log:
        CONCOCT / "run" / "{assembly_id}.log",
    conda:
        "concoct.yml"
    params:
        basename=compose_basename_for_concoct_run_one,
    threads: 24
    resources:
        runtime=24 * 60,
        mem_mb=double_ram(4),
    retries: 5
    shell:
        """
        concoct \
            --threads {threads} \
            --composition_file {input.assembly_10k} \
            --coverage_file {input.coverage} \
            --basename {params.basename} \
        2>> {log} 1>&2
        """


rule _assemble__concoct__merge_cutup_clustering:
    """Run concoct merge_cutup_clustering.py on one assembly."""
    input:
        clustering=CONCOCT / "run" / "{assembly_id}_clustering_gt1000.csv",
    output:
        clustering_merged=CONCOCT / "merge" / "{assembly_id}.csv",
    log:
        CONCOCT / "merge" / "{assembly_id}.log",
    conda:
        "concoct.yml"
    shell:
        """
        merge_cutup_clustering.py \
            {input.clustering} \
        > {output.clustering_merged} \
        2>> {log}
        """


rule _assemble__concoct__extract_fasta_bins:
    """Run concoct extract_fasta_bins.py on one assembly."""
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
        clustering_merged=CONCOCT / "merge" / "{assembly_id}.csv",
    output:
        bins=directory(CONCOCT / "fasta_bins" / "{assembly_id}/"),
    log:
        CONCOCT / "fasta_bins/{assembly_id}.log",
    conda:
        "concoct.yml"
    resources:
        mem_mb=4 * 1024,
    shell:
        """
        mkdir --parents {output.bins}
        extract_fasta_bins.py \
            {input.assembly} \
            {input.clustering_merged} \
            --output_path {output.bins} \
        2> {log} 1>&2
        """


rule assemble__concoct:
    """Run concoct on all assemblies."""
    input:
        [CONCOCT / "fasta_bins" / assembly_id for assembly_id in ASSEMBLIES],
