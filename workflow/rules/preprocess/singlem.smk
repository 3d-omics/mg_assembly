rule preprocess__singlem__fastq:
    input:
        cram=get_host_clean_cram,
    output:
        forward_=temp(SINGLEM / "fastq" / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(SINGLEM / "fastq" / "{sample_id}.{library_id}_2.fq.gz"),
    log:
        SINGLEM / "fastq" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( samtools view \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools collate \
            -O \
            -u \
            -f \
            --threads {threads} \
            - \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 1 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule preprocess__singlem__pipe:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=SINGLEM / "fastq" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=SINGLEM / "fastq" / "{sample_id}.{library_id}_2.fq.gz",
        data=features["databases"]["singlem"],
    output:
        archive_otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "pipe" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --assignment-threads {threads} \
        2> {log} 1>&2
        """


rule preprocess__singlem__condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / "pipe" / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        data=features["databases"]["singlem"],
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule preprocess__singlem__microbial_fraction:
    """Run singlem microbial_fraction over one sample"""
    input:
        forward_=SINGLEM / "fastq" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=SINGLEM / "fastq" / "{sample_id}.{library_id}_2.fq.gz",
        data=features["databases"]["singlem"],
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv",
    output:
        microbial_fraction=SINGLEM
        / "microbial_fraction"
        / "{sample_id}.{library_id}.tsv",
    log:
        SINGLEM / "microbial_fraction" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    shell:
        """
        singlem microbial_fraction \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --input-profile {input.condense} \
            --output-tsv {output.microbial_fraction} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule preprocess__singlem__aggregate_microbial_fraction:
    """Aggregate all the microbial_fraction files into one tsv"""
    input:
        tsvs=[
            SINGLEM / "microbial_fraction" / f"{sample_id}.{library_id}.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=SINGLEM / "microbial_fraction.tsv",
    log:
        SINGLEM / "microbial_fraction.log",
    conda:
        "__environment__.yml"
    shell:
        """
        ( csvstack \
            --tabs \
            {input.tsvs} \
        | csvformat \
            --out-tabs \
        > {output.tsv} \
        ) 2> {log}
        """


rule preprocess__singlem:
    """Run all stats singlem rules"""
    input:
        rules.preprocess__singlem__condense.output,
        rules.preprocess__singlem__aggregate_microbial_fraction.output,
