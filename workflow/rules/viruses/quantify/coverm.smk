rule _viral__quantify__coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        cram=VBOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=VBOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=MMSEQS / "rep_seq.fasta",
        fai=MMSEQS / "rep_seq.fasta.fai",
    output:
        tsv=VCOVERM / "genome" / "{method}" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    log:
        VCOVERM / "genome" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm genome \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output.tsv} \
        ) 2> {log}
        """


rule _viral__quantify__coverm__genome_aggregate:
    """Run coverm genome and a single method"""
    input:
        get_tsvs_for_dereplicate_vcoverm_genome,
    output:
        tsv=VCOVERM / "genome.{method}.tsv",
    log:
        VCOVERM / "genome.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: VCOVERM / "genome" / w.method,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule viral__quantify__coverm__genome:
    """Run coverm genome and all methods"""
    input:
        [
            VCOVERM / f"genome.{method}.tsv"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
        ],


# coverm contig ----
rule _viral__quantify__coverm__contig:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        cram=VBOWTIE2 / "{sample_id}.{library_id}.cram",
        crai=VBOWTIE2 / "{sample_id}.{library_id}.cram.crai",
        reference=MMSEQS / "rep_seq.fasta",
        fai=MMSEQS / "rep_seq.fasta.fai",
    output:
        tsv=VCOVERM / "contig" / "{method}" / "{sample_id}.{library_id}.tsv",
    conda:
        "__environment__.yml"
    log:
        VCOVERM / "contig" / "{method}" / "{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        ( samtools view \
            --exclude-flags 4 \
            --reference {input.reference} \
            --fast \
            {input.cram} \
        | coverm contig \
            --bam-files /dev/stdin \
            --methods {params.method} \
            --proper-pairs-only \
        > {output.tsv} \
        ) 2> {log}
        """


rule _viral__quantify__coverm__contig_aggregate:
    """Run coverm contig and a single method"""
    input:
        get_tsvs_for_dereplicate_vcoverm_contig,
    output:
        tsv=VCOVERM / "contig.{method}.tsv",
    log:
        VCOVERM / "contig.{method}.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=lambda w: VCOVERM / "contig" / w.method,
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        Rscript --vanilla workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule viral__quantify__coverm__contig:
    """Run coverm contig and all methods"""
    input:
        [
            VCOVERM / f"contig.{method}.tsv"
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule viral__quantify__coverm:
    input:
        rules.viral__quantify__coverm__contig.input,
        rules.viral__quantify__coverm__genome.input,
