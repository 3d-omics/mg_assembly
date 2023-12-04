def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    return [
        FASTP / f"{sample_id}.{library_id}_1.fq.gz"
        if len(HOST_NAMES) == 0
        else PRE_BOWTIE2 / f"non{LAST_HOST}" / f"{sample_id}.{library_id}_1.fq.gz"
        for sample_id, library_id in samples_in_assembly
    ]


def get_reverses_from_assembly_id(wildcards):
    """Get the reverse files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    return [
        FASTP / f"{sample_id}.{library_id}_2.fq.gz"
        if len(HOST_NAMES) == 0
        else PRE_BOWTIE2 / f"non{LAST_HOST}" / f"{sample_id}.{library_id}_1.fq.gz"
        for sample_id, library_id in samples_in_assembly
    ]


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = [str(forward_) for forward_ in get_forwards_from_assembly_id(wildcards)]
    return ",".join(forwards)


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = [str(reverse_) for reverse_ in get_reverses_from_assembly_id(wildcards)]
    return ",".join(reverses)


# def get_number_of_libraries_in_assembly(wildcards):
#     """Get the number of libraries in an assembly"""
#     assembly_id = wildcards.assembly_id
#     number_of_libraries = len(samples[samples.assembly_id == assembly_id].index)
#     return number_of_libraries


# def get_crams_to_merge_assembly(wildcards):
#     """Get all the cram files involver in a concrete assembly"""
#     assembly_id = wildcards.assembly_id
#     samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
#     cram_files = [
#         ASSEMBLY_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram"
#         for sample_id, library_id in samples_in_assembly
#     ]
#     return cram_files


# def get_tsvs_for_assembly_coverm_genome(wildcards):
#     """Get all the concrete coverm genome tsv files for a concrete assembly"""
#     tsv_files = [
#         ASSEMBLE_COVERM
#         / "genome"
#         / wildcards.method
#         / f"{assembly_id}.{sample_id}.{library_id}.tsv"
#         for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
#     ]
#     return tsv_files


# def get_tsvs_for_assembly_coverm_contig(wildcards):
#     """Get all the concrete coverm contig tsv files for a concrete assembly"""
#     tsv_files = [
#         ASSEMBLE_COVERM
#         / "contig"
#         / wildcards.method
#         / f"{assembly_id}.{sample_id}.{library_id}.tsv"
#         for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
#     ]
#     return tsv_files


def compose_out_dir_for_assemble_megahit_one(wildcards):
    """Compose output folder"""
    return MEGAHIT / wildcards.assembly_id


# def compose_input_dir_for_assemble_coverm_aggregate_contig(wildcards):
#     """Compose the input dir for coverme contig"""
#     return ASSEMBLE_COVERM / "contig" / wildcards.method


# def compose_input_dir_for_assemble_coverm_aggregate_genome(wildcards):
#     """Compose the input dir for coverm genome"""
#     return ASSEMBLE_COVERM / "genome" / wildcards.method


def get_bams_from_assembly_id(wildcards):
    """Given an assembly_id, get all the bam files for that assembly."""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = [
        ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam"
        for sample_id, library_id in samples_in_assembly
    ]
    return bam_files


def get_bais_from_assembly_id(wildcards):
    """Given an assembly_id, get all the bai files for that assembly."""
    bams = get_bams_from_assembly_id(wildcards)
    return [f"{bam}.bai" for bam in bams]


def compose_out_prefix_for_bin_magscot_run_one(wildcards):
    """Compose the output folder for magscot"""
    return MAGSCOT / wildcards.assembly_id / "magscot"


def compose_basename_for_concoct_run_one(wildcards):
    """Compose the basename for the concoct run"""
    return CONCOCT / "run" / wildcards.assembly_id


def compose_out_prefix_for_maxbin2_run_one(wildcards):
    """Compose the output folder for maxbin2"""
    return MAXBIN2 / "bins" / wildcards.assembly_id


def compose_bins_prefix_for_metabat2_run_one(wildcards):
    """Compose the output folder for metabat2"""
    return METABAT2 / "bins" / wildcards.assembly_id / wildcards.assembly_id
