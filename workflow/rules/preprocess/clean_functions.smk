# functions to choose between fastp or last host
def get_final_fastq(wildcards):
    """Get the final fastq file for the host mapping"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return [
            PRE_FASTP / f"{sample_id}.{library_id}_1.fq.gz",
            PRE_FASTP / f"{sample_id}.{library_id}_2.fq.gz",
        ]
    last_host = HOST_NAMES[-1]
    return [
        PRE_BOWTIE2 / f"{last_host}" / f"{sample_id}.{library_id}_u1.fq.gz",
        PRE_BOWTIE2 / f"{last_host}" / f"{sample_id}.{library_id}_u2.fq.gz",
    ]


def get_final_fastq_forward(wildcards):
    return get_final_fastq(wildcards)[0]


def get_final_fastq_reverse(wildcards):
    return get_final_fastq(wildcards)[1]
