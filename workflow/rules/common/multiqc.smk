rule multiqc:
    """Aggregate QC reports with MultiQC"""
    input:
        "results",
    output:
        "report.html",
        "report_data.zip",
    log:
        "report.log",
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=6 * 60,
    params:
        extra="--title report --dirs --fullnames --fn_as_s_name --force",
    wrapper:
        "v9.4.0/bio/multiqc"
