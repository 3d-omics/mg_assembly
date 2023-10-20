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
        cram_files.append(BIN_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.cram")
    return cram_files


def get_tsvs_for_binning_coverm_genome(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            BIN_COVERM / "genome" / f"{assembly_id}.{sample_id}.{library_id}.tsv"
        )
    return tsv_files


def get_tsvs_for_binning_coverm_contig(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            BIN_COVERM / "contig" / f"{assembly_id}.{sample_id}.{library_id}.tsv"
        )
    return tsv_files


# Magscot ----
def compose_out_prefix_for_metabin_magscot_run_one(wildcards):
    return MAGSCOT / wildcards.assembly_id / "magscot"


# Concoct ----
def get_bams_for_concoct_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files


def get_bais_for_concoct_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam.bai",
        )
    return bam_files


def get_basename_for_concoct_run_one(wildcards):
    return CONCOCT / "run" / wildcards.assembly_id


# MaxBin2 ----
def get_bams_for_maxbin2(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files


def compose_out_prefix_for_maxbin2_run_one(wildcards):
    return MAXBIN2 / "bins" / wildcards.assembly_id


# MetaBat2 ---
def get_bams_for_metabat2(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files


def compose_bins_prefix_for_metaba2_run_one(wildcards):
    return METABAT2 / "bins" / wildcards.assembly_id / wildcards.assembly_id
