def get_bowtie2_host_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_hosts reports"""
    sample = wildcards.sample
    library = wildcards.library
    return [
        BOWTIE2_HOSTS / f"{host}/{sample}.{library}.{report}"
        for host in ["human", "chicken"]
        for report in BAM_REPORTS
    ]


def get_bowtie2_mags_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_mags reports"""
    sample = wildcards.sample
    library = wildcards.library
    return [
        BOWTIE2_MAGS / f"{mag_catalogue}/{sample}.{library}.{report}"
        for mag_catalogue in MAG_CATALOGUES
        for report in BAM_REPORTS
    ]


def get_kraken2_for_library_reports(wildcards):
    """Compose the paths for the kraken2 reports"""
    sample = wildcards.sample
    library = wildcards.library
    return [
        f"{KRAKEN2}/{kraken2_db}/{sample}.{library}.report"
        for kraken2_db in KRAKEN2_DBS
    ]
