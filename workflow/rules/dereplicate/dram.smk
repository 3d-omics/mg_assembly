rule dereplicate_dram_setup_db:
    """Set up the DRAM database."""
    input:
        features["dram_database"],
    output:
        touch(DREP_DRAM / "dram_db_setup.done"),
    log:
        DREP_DRAM / "dram_db_setup.log",
    conda:
        "dram.yml"
    shell:
        "python workflow/scripts/dram_setup.py {input} 2> {log}"


rule dereplicate_dram_annotate:
    """Annotate dereplicated genomes with DRAM."""
    input:
        drep_folder=DREP / "dereplicated_genomes",
        mock_db=DREP_DRAM / "dram_db_setup.done",
    output:
        outdir=DREP_DRAM / "annotate" / "dereplicated_genomes",
        annotations=DREP_DRAM / "annotate" / "dereplicated_genomes" / "annotations.tsv",
        trnas=touch(DREP_DRAM / "annotate" / "dereplicated_genomes" / "trnas.tsv"),
        rrnas=touch(DREP_DRAM / "annotate" / "dereplicated_genomes" / "rrnas.tsv"),
    log:
        DREP_DRAM / "annotate" / "dereplicated_genomes.log",
    conda:
        "dram.yml"
    params:
        min_contig_size=1500,
    threads: 24
    resources:
        mem_mb=double_ram(params["dereplicate"]["dram"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        rm -rfv {output.outdir} 2> {log} 1>&2

        DRAM.py annotate \
            --input_fasta {input.drep_folder} \
            --output_dir {output.outdir} \
            --threads {threads} \
            --rrna_path {output.rrnas} \
            --trna_path {output.trnas} \
            --min_contig_size {params.min_contig_size} \
        2>> {log} 1>&2
        """


rule dereplicate_dram_distill:
    """Distill DRAM annotations."""
    input:
        indir=DREP_DRAM / "annotate" / "dereplicated_genomes",
        annotations=DREP_DRAM / "annotate" / "dereplicated_genomes" / "annotations.tsv",
        trnas=DREP_DRAM / "annotate" / "dereplicated_genomes" / "trnas.tsv",
        rrnas=DREP_DRAM / "annotate" / " dereplicated_genomes" / "rrnas.tsv",
        mock_db=DREP_DRAM / "dram_db_setup.done",
    output:
        outdir=DREP_DRAM / "distill" / "dereplicated_genomes",
    log:
        DREP_DRAM / "distill" / "dereplicated_genomes.log",
    conda:
        "dram.yml"
    resources:
        mem_mb=double_ram(params["dereplicate"]["dram"]["memory_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        DRAM.py distill \
            --input_dir {input.indir} \
            --output_dir {output.outdir} \
            --annotations {input.annotations} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
        2> {log} 1>&2
        """


rule dereplicate_dram:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.dereplicate_dram_distill.output,


localrules:
    dereplicate_dram_setup_db,
