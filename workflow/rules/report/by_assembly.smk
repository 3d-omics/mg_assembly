rule report_assembly_one:
    """Create a report for a single assembly

    To do so, collect all the files and all the libraries related to that assembly
    """
    input:
        get_stats_files_from_assembly_id,
    output:
        report=REPORT_ASSEMBLY / "{assembly_id}.html",
    log:
        REPORT_ASSEMBLY / "{assembly_id}.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_ASSEMBLY,
        title="{assembly_id}",
    shell:
        """
        multiqc \
            --title {params.title} \
            --force \
            --filename {params.title} \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report_assembly:
    input:
        [REPORT_ASSEMBLY / f"{assembly_id}.html" for assembly_id in ASSEMBLIES],
