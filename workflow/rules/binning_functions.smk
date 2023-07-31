# def get_assembly_bams_from_assembly_id(wildcards):
#     assembly_id = wildcards.assembly_id
#     samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
#     bam_files = []
#     for sample_id, library_id in samples_in_assembly:
#         bam_files.append(
#             BOWTIE2_ASSEMBLY / f"{assembly_id}.{sample_id}.{library_id}.bam"
#         )
#     return bam_files


def get_number_of_libraries_in_binning(wildcards):
    """Get the number of libraries in an assembly"""
    assembly_id = wildcards.assembly_id
    number_of_libraries = len(samples[samples.assembly_id == assembly_id].index)
    return number_of_libraries


def get_crams_to_merge_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    cram_files = []
    for sample_id, library_id in samples_in_assembly:
        cram_files.append(
            BOWTIE2_BINNING / f"{assembly_id}.{sample_id}.{library_id}.cram"
        )
    return cram_files


def get_tsvs_for_binning_coverm_genome(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            COVERM_BINNING / f"genome/{assembly_id}.{sample_id}.{library_id}.tsv"
        )
    return tsv_files


def get_tsvs_for_binning_coverm_contig(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            COVERM_BINNING / f"contig/{assembly_id}.{sample_id}.{library_id}.tsv"
        )
    return tsv_files
