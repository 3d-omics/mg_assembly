rule pre_singlem_pipe_one:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=get_final_forward_from_pre,
        reverse_=get_final_reverse_from_pre,
        data=features["singlem_database"],
    output:
        archive_otu_table=SINGLEM / "{sample_id}.{library_id}.archive.json",
        otu_table=SINGLEM / "{sample_id}.{library_id}.otu_table.tsv",
        condense=SINGLEM / "{sample_id}.{library_id}.condense.tsv",
    log:
        SINGLEM / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
    threads: 24
    resources:
        runtime=4 * 60,
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --threads {threads} \
            --assignment-threads {threads} \
        2> {log} 1>&2 || true
        """


rule pre_singlem:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / f"{sample_id}.{library_id}.archive.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=features["singlem_database"],
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    conda:
        "_env.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.database} \
        2> {log} 1>&2
        """
