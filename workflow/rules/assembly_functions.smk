def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    assembly_id = wildcards.assembly_id
    forward_filenames = samples[
        samples.assembly_id == assembly_id
    ].forward_filename.values.tolist()
    return forward_filenames


def get_reverses_from_assembly_id(wildcards):
    """Get the reverse files for megahit"""
    assembly_id = wildcards.assembly_id
    reverse_filenames = samples[
        samples.assembly_id == assembly_id
    ].reverse_filename.values.tolist()
    return reverse_filenames


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = get_forwards_from_assembly_id(wildcards)
    return ",".join(forwards)


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = get_reverses_from_assembly_id(wildcards)
    return ",".join(reverses)


def get_number_of_libraries_in_assembly(wildcards):
    """Get the number of libraries in an assembly"""
    assembly_id = wildcards.assembly_id
    number_of_libraries = len(samples[samples.assembly_id == assembly_id].index)
    return number_of_libraries


def get_crams_to_merge(wildcards):
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


def compose_metawrap_working_folder(wildcards):
    assembly_id = wildcards.assembly_id
    completeness = params["assembly"]["metawrap_bin_refinement"]["completeness"]
    contamination = params["assembly"]["metawrap_bin_refinement"]["contamination"]
    folder = (
        METAWRAP_REFINEMENT
        / f"{assembly_id}/metawrap_{completeness}_{contamination}_bins"
    )
    return folder
