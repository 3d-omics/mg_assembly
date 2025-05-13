# generic functions to infer sample and library ids
def compose_rg_id(wildcards):
    """Compose the read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose the read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.sample_id}"


# functions to choose between fastp or already mapped files
def get_fastq_for_host_mapping(wildcards):
    """Get the input cram file for host mapping"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    host = wildcards.host
    host_index = HOST_NAMES.index(host)
    if host_index == 0:
        return [
            PRE_FASTP / f"{sample_id}.{library_id}_1.fq.gz",
            PRE_FASTP / f"{sample_id}.{library_id}_2.fq.gz",
        ]
    previous_host = HOST_NAMES[host_index - 1]
    return [
        PRE_BOWTIE2 / f"{previous_host}" / f"{sample_id}.{library_id}_u1.fq.gz",
        PRE_BOWTIE2 / f"{previous_host}" / f"{sample_id}.{library_id}_u2.fq.gz",
    ]


def get_fastq_for_host_mapping_forward(wildcards):
    return get_fastq_for_host_mapping(wildcards)[0]


def get_fastq_for_host_mapping_reverse(wildcards):
    return get_fastq_for_host_mapping(wildcards)[1]
