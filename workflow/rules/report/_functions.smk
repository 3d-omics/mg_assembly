# pre -----
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


def get_sample_library_from_assembly(my_assembly_id):
    """Get all the samples involved in a given assembly"""
    SAMPLE_LIBRARY = [
        [sample_id, library_id]
        for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        if assembly_id == my_assembly_id
    ]
    return SAMPLE_LIBRARY


def _path_to_string(path_list):
    """Convert a list of paths to a list of strings"""
    return [str(x) for x in path_list]


def get_stats_files_from_assembly_id(wildcards):
    """Compose the paths for the reads fastqc reports"""
    sample_library = [
        [sample_id, library_id]
        for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        if assembly_id == wildcards.assembly_id
    ]

    # reads ----
    reads_fastqc = [
        READS / f"{sample}.{library}_{end}_fastqc.zip"
        for sample, library in sample_library
        for end in ["1", "2"]
    ]

    # preprocessing ----
    pre_fastp_fastqc = [
        FASTP / f"{sample}.{library}_{end}_fastqc.zip"
        for sample, library in sample_library
        for end in ["1", "2", "u1", "u2"]
    ]
    pre_fastp = [
        FASTP / f"{sample_id}.{library_id}_fastp.html"
        for sample_id, library_id in sample_library
    ]
    pre_samtools = (
        [
            PRE_BOWTIE2 / genome / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in sample_library
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )
    pre_nonhost_fastqc = (
        [
            PRE_BOWTIE2 / f"non{genome}/{sample_id}.{library_id}_{end}_fastqc.zip"
            for sample_id, library_id in sample_library
            for end in ["1", "2"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )
    pre_kraken2 = [
        KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
        for sample_id, library_id in sample_library
        for kraken_db in KRAKEN2_DBS
    ]

    # assembly ----
    assemble_quast = ASSEMBLE_QUAST / f"{wildcards.assembly_id}"
    assemble_samtools = [
        ASSEMBLE_BOWTIE2
        / f"{wildcards.assembly_id}.{sample_id}.{library_id}.{extension}"
        for sample_id, library_id in sample_library
        for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
    ]

    all_files = (
        _path_to_string(reads_fastqc)
        + _path_to_string(pre_fastp_fastqc)
        + _path_to_string(pre_fastp)
        + _path_to_string(pre_samtools)
        + _path_to_string(pre_nonhost_fastqc)
        + _path_to_string(pre_kraken2)
        + _path_to_string(assemble_samtools)
        + [str(assemble_quast)]
    )

    return all_files


def get_stats_files_from_sample_and_library_ids(wildcards):
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id

    reads_fastqc = [
        READS / f"{sample_id}.{library_id}_{end}_fastqc.zip" for end in ["1", "2"]
    ]

    pre_fastp_fastqc = [
        FASTP / f"{sample_id}.{library_id}_{end}_fastqc.zip"
        for end in ["1", "2", "u1", "u2"]
    ]

    pre_fastp_html = FASTP / f"{sample_id}.{library_id}_fastp.html"

    pre_bowtie2 = (
        [
            PRE_BOWTIE2 / genome / f"{sample_id}.{library_id}.{extension}"
            for extension in ["stats.txt", "flagstats.txt", "idxstats.tsv"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )

    pre_nonhost_fastqc = (
        [
            PRE_BOWTIE2 / f"non{genome}" / f"{sample_id}.{library_id}_{end}_fastqc.zip"
            for end in ["1", "2"]
            for genome in HOST_NAMES
        ]
        if len(HOST_NAMES) > 0
        else []
    )

    pre_kraken2 = [
        KRAKEN2 / kraken_db / f"{sample_id}.{library_id}.report"
        for kraken_db in KRAKEN2_DBS
    ]

    quantify_bowtie2 = [
        QUANT_BOWTIE2 / f"{sample_id}.{library_id}.{extension}"
        for extension in ["stats.txt", "flagstats.txt"]
    ]

    all_files = (
        _path_to_string(reads_fastqc)
        + _path_to_string(pre_fastp_fastqc)
        + [pre_fastp_html]
        + _path_to_string(pre_bowtie2)
        + _path_to_string(pre_nonhost_fastqc)
        + _path_to_string(pre_kraken2)
        + _path_to_string(quantify_bowtie2)
    )

    return all_files
