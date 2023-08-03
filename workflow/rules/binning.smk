include: "binning/vamb.smk"
include: "binning/concoct.smk"
include: "binning/metabat2.smk"
include: "binning/maxbin2.smk"
include: "binning/metawrap.smk"
include: "binning/magscot.smk"


rule binning_index_one:
    input:
        bins=MAGSCOT / "{assembly_id}.fa",
    output:
        mock=touch(METABIN_INDEX / "{assembly_id}"),
    log:
        METABIN_INDEX / "{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    params:
        extra=params["assembly"]["bowtie2-build"]["extra"],
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.bins} \
            {output.mock} \
        2> {log} 1>&2
        """


rule binning_bowtie2_one:
    input:
        mock=METABIN_INDEX / "{assembly_id}",
        forward_=NONHOST / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=NONHOST / "{sample_id}.{library_id}_2.fq.gz",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        cram=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.cram",
    log:
        METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 24
    params:
        extra=params["binning"]["bowtie2"]["extra"],
        samtools_mem=params["binning"]["samtools"]["mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    resources:
        mem_mb=32 * 1024,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule binning_bowtie2:
    input:
        [
            METABIN_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule binning_cram_to_bam_one:
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
        "../envs/binning.yml"
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


rule binning_coverm_genome_one:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        tsv=METABIN_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/binning.yml"
    log:
        METABIN_COVERM / "genome/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["assembly"]["coverm"]["genome"]["methods"],
        min_covered_fraction=params["assembly"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["assembly"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} 2> {log}
        """


rule binning_coverm_genome:
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


rule binning_coverm_contig_one:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=METABIN_BOWTIE2 / "{assembly_id}.{sample_id}.{library_id}.bam",
        reference=MAGSCOT / "{assembly_id}.fa",
    output:
        tsv=METABIN_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.tsv",
    conda:
        "../envs/binning.yml"
    log:
        METABIN_COVERM / "contig/{assembly_id}.{sample_id}.{library_id}.log",
    params:
        methods=params["binning"]["coverm"]["contig"]["methods"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output.tsv} 2> {log}
        """


rule binning_coverm_contig:
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
        "../envs/binning.yml"
    params:
        input_dir=METABIN_COVERM / f"contig",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule binning_quast_one:
    """Run quast over one assembly group"""
    input:
        MAGSCOT / "{assembly_id}.fa",
    output:
        directory(METABIN_QUAST / "{assembly_id}"),
    log:
        METABIN_QUAST / "{assembly_id}.log",
    conda:
        "../envs/binning.yml"
    threads: 4
    params:
        extra=params["assembly"]["quast"]["extra"],
    shell:
        """
        quast \
            --output-dir {output} \
            --threads {threads} \
            {params.extra} \
            {input} \
        2> {log} 1>&2
        """


rule binning_quast_all:
    """Run quast over all assembly groups"""
    input:
        [METABIN_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule binning_gtdbtk_one:
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


rule binning_gtdbtk:
    input:
        [METABIN_GTDBTK / f"{assembly_id}" for assembly_id in ASSEMBLIES],


rule binning_dram_annotate_one:
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
        "../envs/binning.yml"
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


rule binning_dram_distill_one:
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
        "../envs/binning.yml"
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


rule binning_dram:
    input:
        [METABIN_DRAM / f"distill/{assembly_id}" for assembly_id in ASSEMBLIES],


rule binning:
    input:
        rules.binning_coverm_contig.output,
        rules.binning_coverm_genome.output,
        rules.binning_quast_all.input,
