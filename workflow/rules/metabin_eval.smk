rule metabin_eval_cram_to_bam_one:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
        crai=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram.crai",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        bam=temp(METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam"),
    log:
        METABIN_BOWTIE2 / "{assembly_id},{sample_id}.{library_id}.bam.log",
    conda:
        "../envs/metabin.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --threads {threads} \
            --reference {input.reference} \
            --output {output.bam} \
            --fast \
            {input.cram} \
        2> {log}
        """


rule metabin_eval_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        tsv=METABIN_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/metabin.yml"
    log:
        METABIN_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["metabin"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["metabin"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["metabin"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule metabin_eval_coverm_genome:
    input:
        tsvs=[
            METABIN_COVERM / f"genome/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=METABIN_COVERM / "genome.tsv",
    log:
        METABIN_COVERM / "genome.log",
    conda:
        "../envs/assembly.yml"
    params:
        input_dir=METABIN_COVERM / "genome/",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule metabin_eval_coverm_contig_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        tsv=METABIN_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/metabin.yml"
    log:
        METABIN_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["metabin"]["coverm"]["contig"]["methods"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule metabin_eval_coverm_contig:
    input:
        tsvs=[
            METABIN_COVERM / f"contig/{assembly_id}.{sample_id}.{library_id}.tsv"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
    output:
        tsv=METABIN_COVERM / "contig.tsv",
    log:
        METABIN_COVERM / "contig.log",
    conda:
        "../envs/metabin.yml"
    params:
        input_dir=METABIN_COVERM / f"contig",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule metabin_eval_quast_one:
    """Run quast over one assembly group"""
    input:
        MAGSCOT / "{assembly_id}.fa",
    output:
        directory(METABIN_QUAST / "{assembly_id}"),
    log:
        METABIN_QUAST / "{assembly_id}.log",
    conda:
        "../envs/metabin.yml"
    threads: 4
    params:
        extra=params["metabin"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule metabin_eval_quast_all:
    """Run quast over all assembly groups"""
    input:
        [METABIN_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule metabin_eval_gtdbtk_one:
    input:
        bin_folder=MAGSCOT / "{assembly_id}/bins",
        database=features["gtdbtk_database"],
    output:
        outdir=METABIN_GTDBTK / "{assembly_id}",
    threads: 16
    resources:
        mem_mb=150 * 1024,
    shell:
        """
        GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {params.bins} \
            --extension gz \
            --out_dir {output.outdir} \
            --cpus {threads} \
            --skip_ani_screen \
        2> {log} 1>&2
        """


rule metabin_eval_gtdbtk:
    input:
        [METABIN_GTDBTK / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule metabin_eval_dram_annotate_one:
    input:
        bin_folder=MAGSCOT / "{assembly_id}/bins/",
        dram_database=features["dram_database"],
    output:
        outdir=directory(METABIN_DRAM / "annotate/{assembly_id}"),
        annotations=METABIN_DRAM / "annotate/{assembly_id}/annotations.tsv",
        trnas=touch(METABIN_DRAM / "annotate/{assembly_id}/trnas.tsv"),
        rrnas=touch(METABIN_DRAM / "annotate/{assembly_id}/rrnas.tsv"),
    log:
        METABIN_DRAM / "annotate/{assembly_id}.log",
    conda:
        "../envs/metabin.yml"
    params:
        min_contig_size=1500,
    shell:
        """
        DRAM.py annotate \
            --input_fasta {input.bin_fa} \
            --output_dir {params.outdir} \
            --threads {threads} \
            --min_contig_size {params.min_contig_size}
        """


rule metabin_eval_dram_distill_one:
    input:
        outdir=METABIN_DRAM / "annotate/{assembly_id}",
        annotations=METABIN_DRAM / "annotate/{assembly_id}/annotations.tsv",
        trnas=METABIN_DRAM / "annotate/{assembly_id}/trnas.tsv",
        rrnas=METABIN_DRAM / "annotate/{assembly_id}/rrnas.tsv",
    output:
        outdir=directory(METABIN_DRAM / "distill/{assembly_id}"),
    log:
        METABIN_DRAM / "distill/{assembly_id}.log",
    conda:
        "../envs/metabin.yml"
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --output_file {output.outdir} \
            --threads {threads} \
            --rrna_path {input.rrnas} \
            --trna_path {input.trnas} \
        2> {log} 1>&2
        """


rule metabin_eval_dram:
    input:
        [METABIN_DRAM / f"distill/{assembly_id}" for assembly_id in ASSEMBLIES],


rule metabin_eval:
    input:
        rules.metabin_eval_coverm_contig.output,
        rules.metabin_eval_coverm_genome.output,
        rules.metabin_eval_quast_all.input,
        rules.metabin_eval_dram.input,
        rules.metabin_eval_gtdbtk.input,
