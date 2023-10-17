def get_tsvs_for_dereplicate_coverm_genome(wildcards):
    method = wildcards.method
    tsv_files = [
        DEREPLICATE_COVERM / f"genome/{method}/{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_dereplicate_coverm_contig(wildcards):
    method = wildcards.method
    tsv_files = [
        DEREPLICATE_COVERM / f"contig/{method}/{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def double_ram_for_dereplicate_drep(wildcards, attempt):
    initial_memory = params["dereplicate"]["drep"]["memory_gb"]
    return initial_memory * 2**attempt * 1024


def double_ram_for_dereplicate_bowtie2_build(wildcards, attempt):
    initial_memory = params["dereplicate"]["bowtie2-build"]["memory_gb"]
    return initial_memory * 2**attempt * 1024


def double_ram_for_dereplicate_bowtie2_map(wildcards, attempt):
    initial_memory = params["dereplicate"]["bowtie2"]["memory_gb"]
    return initial_memory * 2**attempt * 1024


def double_ram_for_dereplicate_eval_gtdbtk(wildcards, attempt):
    initial_memory = params["dereplicate"]["gtdbtk"]["memory_gb"]
    return initial_memory * 2**attempt * 1024


def double_ram_for_dereplicate_eval_dram(wildcards, attempt):
    initial_memory = params["dereplicate"]["dram"]["memory_gb"]
    return initial_memory * 2**attempt * 1024
