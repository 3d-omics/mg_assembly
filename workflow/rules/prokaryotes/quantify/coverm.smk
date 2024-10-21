rule prokaryotes__quantify__coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.bam",
        bai=QUANT_BOWTIE2
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.bam.bai",
        reference=PROK_ANN / "drep.{secondary_ani}.fa.gz",
        fai=PROK_ANN / "drep.{secondary_ani}" / "dereplicated_genomes.fa.gz.fai",
    output:
        tsv=COVERM
        / "genome"
        / "{method}"
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        COVERM
        / "genome"
        / "{method}"
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        ( coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        | gzip --best \
        > {output.tsv} \
        ) 2> {log}
        """


rule prokaryotes__quantify__coverm__genome__aggregate:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_genome,
    output:
        tsv=COVERM / "genome.{method}.{secondary_ani}.tsv.gz",
    log:
        COVERM / "genome.{method}.{secondary_ani}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: COVERM / "genome" / w.method / w.secondary_ani,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule prokaryotes__quantify__coverm__genome__all:
    """Run coverm genome and all methods"""
    input:
        [
            COVERM / f"genome.{method}.{secondary_ani}.tsv.gz"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
            for secondary_ani in SECONDARY_ANIS
        ],


# coverm contig ----
rule prokaryotes__quantify__coverm__contig:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=QUANT_BOWTIE2 / "drep.{secondary_ani}" / "{sample_id}.{library_id}.bam",
        bai=QUANT_BOWTIE2
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.bam.bai",
        reference=PROK_ANN / "drep.{secondary_ani}.fa.gz",
        fai=PROK_ANN / "drep.{secondary_ani}.fa.gz.fai",
    output:
        tsv=COVERM
        / "contig"
        / "{method}"
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.tsv.gz",
    conda:
        "__environment__.yml"
    log:
        COVERM
        / "contig"
        / "{method}"
        / "drep.{secondary_ani}"
        / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        ( coverm contig \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --proper-pairs-only \
        | gzip --best \
        > {output.tsv} \
        ) 2> {log}
        """


rule prokaryotes__quantify__coverm__contig__aggregate:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_coverm_contig,
    output:
        tsv=COVERM / "contig.{method}.{secondary_ani}.tsv.gz",
    log:
        COVERM / "contig.{method}.{secondary_ani}log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: COVERM / "contig" / w.method / w.secondary_ani,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule prokaryotes__quantify__coverm__contig__all:
    """Run coverm contig and all methods"""
    input:
        [
            COVERM / f"contig.{method}.{secondary_ani}.tsv.gz"
            for method in params["quantify"]["coverm"]["contig"]["methods"]
            for secondary_ani in SECONDARY_ANIS
        ],


rule prokaryotes__quantify__coverm__all:
    input:
        rules.prokaryotes__quantify__coverm__genome__all.input,
        rules.prokaryotes__quantify__coverm__contig__all.input,
