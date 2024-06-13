rule preprocess__singlem__pipe__:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        metapackage=features["databases"]["singlem"],
    output:
        archive_otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.archive.json",  # don't compress
        otu_table=SINGLEM / "pipe" / "{sample_id}.{library_id}.otu_table.tsv.gz",
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv.gz",
    log:
        SINGLEM / "pipe" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=lambda w: SINGLEM / "pipe" / f"{w.sample_id}.{w.library_id}.condense.tsv",
        otu_table=lambda w: SINGLEM
        / "pipe"
        / f"{w.sample_id}.{w.library_id}.otu_table.tsv",
    shell:
        """
        {{
            singlem pipe \
                --forward {input.forward_} \
                --reverse {input.reverse_} \
                --otu-table {output.otu_table} \
                --archive-otu-table {output.archive_otu_table} \
                --taxonomic-profile {output.condense} \
                --metapackage {input.metapackage} \
                --threads {threads} \
                --assignment-threads {threads} \

            gzip \
                --verbose \
                {params.prefix} \
                {params.otu_table}

        }} 2> {log} 1>&2 || true
        """


rule preprocess__singlem__condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / "pipe" / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=features["databases"]["singlem"],
    output:
        condense=SINGLEM / "singlem.tsv.gz",
    log:
        SINGLEM / "singlem.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --taxonomic-profile >(gzip > {output.condense}) \
            --metapackage {input.database} \
            --input-archive-otu-tables {input.archive_otu_tables} \
        2> {log} 1>&2
        """


rule preprocess__singlem__microbial_fraction__:
    """Run singlem microbial_fraction over one sample"""
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        data=features["databases"]["singlem"],
        condense=SINGLEM / "pipe" / "{sample_id}.{library_id}.condense.tsv.gz",
    output:
        microbial_fraction=SINGLEM
        / "microbial_fraction"
        / "{sample_id}.{library_id}.tsv.gz",
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
            --output-tsv >(gzip > {output.microbial_fraction}) \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule preprocess__singlem__aggregate_microbial_fraction:
    """Aggregate all the microbial_fraction files into one tsv"""
    input:
        tsvs=[
            SINGLEM / "microbial_fraction" / f"{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        tsv=SINGLEM / "microbial_fraction.tsv.gz",
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
        | gzip --verbose \
        > {output.tsv} \
        ) 2> {log}
        """


rule preprocess__singlem:
    input:
        rules.preprocess__singlem__condense.output,
        rules.preprocess__singlem__aggregate_microbial_fraction.output,
