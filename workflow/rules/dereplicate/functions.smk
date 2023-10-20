def get_forward_for_dereplicate_bowtie2_one(wildcards):
    if len(HOST_NAMES) > 1:
        return PRE_BOWTIE2 / f"non{LAST_HOST}/{sample_id}.{library_id}_{end}.fq.gz"
    return FASTP / "{sample_id}.{library_id}_1.fq.gz"


def get_reverse_for_dereplicate_bowtie2_one(wildcards):
    if len(HOST_NAMES) > 1:
        return PRE_BOWTIE2 / f"non{LAST_HOST}/{sample_id}.{library_id}_{end}.fq.gz"
    return FASTP / "{sample_id}.{library_id}_2.fq.gz"


def get_tsvs_for_dereplicate_coverm_genome(wildcards):
    method = wildcards.method
    tsv_files = [
        DREP_COVERM / f"genome/{method}/{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_dereplicate_coverm_contig(wildcards):
    method = wildcards.method
    tsv_files = [
        DREP_COVERM / f"contig/{method}/{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files
