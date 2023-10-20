def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    forward_filenames = []
    for sample_id, library_id in samples_in_assembly:
        if len(HOST_NAMES) == 0:
            forward_filenames.append(FASTP / f"{sample_id}.{library_id}_1.fq.gz")
        else:
            forward_filenames.append(
                PRE_BOWTIE2 / f"non{LAST_HOST}/{sample_id}.{library_id}_1.fq.gz"
            )
    return forward_filenames


def get_reverses_from_assembly_id(wildcards):
    """Get the reverse files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    reverse_filenames = []
    for sample_id, library_id in samples_in_assembly:
        if len(HOST_NAMES) == 0:
            reverse_filenames.append(FASTP / f"{sample_id}.{library_id}_2.fq.gz")
        else:
            reverse_filenames.append(
                PRE_BOWTIE2 / f"non{LAST_HOST}/{sample_id}.{library_id}_2.fq.gz"
            )
    return reverse_filenames


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = [str(forward_) for forward_ in get_forwards_from_assembly_id(wildcards)]
    return ",".join(forwards)


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = [str(reverse_) for reverse_ in get_reverses_from_assembly_id(wildcards)]
    return ",".join(reverses)


def get_number_of_libraries_in_assembly(wildcards):
    """Get the number of libraries in an assembly"""
    assembly_id = wildcards.assembly_id
    number_of_libraries = len(samples[samples.assembly_id == assembly_id].index)
    return number_of_libraries


def get_crams_to_merge_assembly(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    cram_files = []
    for sample_id, library_id in samples_in_assembly:
        cram_files.append(
            ASSEMBLY_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
        )
    return cram_files


def get_tsvs_for_assembly_coverm_genome(wildcards):
    method = wildcards.method
    tsv_files = [
        ASSEMBLE_COVERM
        / "genome"
        / method
        / f"{assembly_id}.{sample_id}.{library_id}.tsv"
        for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_assembly_coverm_contig(wildcards):
    method = wildcards.method
    tsv_files = [
        ASSEMBLE_COVERM
        / "contig"
        / method
        / f"assembly_id.{sample_id}.{library_id}.tsv"
        for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
    ]
    return tsv_files


def compose_out_dir_for_assemble_megahit_one(wildcards):
    return MEGAHIT / wildcards.assembly_id


def get_memory_bytes_for_assemble_megahit_one(wildcards, resources):
    # https://github.com/snakemake/snakemake/issues/499
    return resources.mem_mb * 1024**2


def compose_input_dir_for_assemble_eval_coverm_aggregate_contig(wildcards):
    return ASSEMBLE_COVERM / "contig" / wildcards.method


def compose_input_dir_for_assemble_eval_coverm_aggregate_genome(wildcards):
    return ASSEMBLE_COVERM / "genome" / wildcards.method
