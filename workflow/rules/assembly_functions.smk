def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    forward_filenames = []
    for sample_id, library_id in samples_in_assembly:
        forward_filenames.append(NONHOST / f"{sample_id}.{library_id}_1.fq.gz")
    return forward_filenames


def get_reverses_from_assembly_id(wildcards):
    """Get the reverse files for megahit"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    reverse_filenames = []
    for sample_id, library_id in samples_in_assembly:
        reverse_filenames.append(NONHOST / f"{sample_id}.{library_id}_1.fq.gz")
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
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    cram_files = []
    for sample_id, library_id in samples_in_assembly:
        cram_files.append(
            BOWTIE2_ASSEMBLY / f"{assembly_id}/{sample_id}.{library_id}.cram"
        )
    return cram_files


def get_tsvs_for_assembly_coverm_genome(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            COVERM_ASSEMBLY / f"genome/{assembly_id}/{sample_id}.{library_id}.tsv"
        )
    return tsv_files


def get_tsvs_for_assembly_coverm_contig(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    tsv_files = []
    for sample_id, library_id in samples_in_assembly:
        tsv_files.append(
            COVERM_ASSEMBLY / f"contig/{assembly_id}/{sample_id}.{library_id}.tsv"
        )
    return tsv_files
