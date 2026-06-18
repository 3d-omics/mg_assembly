rule viruses__mvp__00_set_up:
    """Validate the MVP metadata and create the MVP working directory skeleton"""
    input:
        VIR_MVP / "metadata.tsv",
    output:
        VIR_MVP / "MVP_00_Summary_Report.txt",
    log:
        VIR_MVP / "mvp_00_set_up.log",
    conda:
        ENVS / "mvp.yml"
    shell:
        """
        mvip MVP_00_set_up_MVP \
            --input {VIR_MVP} \
            --metadata {input} \
            --skip_install_databases \
        2> {log} 1>&2
        """
