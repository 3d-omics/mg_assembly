def compose_metawrap_working_folder(wildcards):
    assembly_id = wildcards.assembly_id
    completeness = params["binning"]["metawrap_bin_refinement"]["completeness"]
    contamination = params["binning"]["metawrap_bin_refinement"]["contamination"]
    folder = (
        METAWRAP_REFINEMENT
        / f"{assembly_id}/metawrap_{completeness}_{contamination}_bins"
    )
    return folder


def get_number_of_libraries_in_binning(wildcards):
    """Get the number of libraries in an assembly"""
    assembly_id = wildcards.assembly_id
    number_of_libraries = len(samples[samples.assembly_id == assembly_id].index)
    return number_of_libraries


def get_crams_to_merge_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    cram_files = []
    for sample_id, library_id in samples_in_assembly:
        cram_files.append(
            BOWTIE2_BINNING / f"{assembly_id}/{sample_id}.{library_id}.cram"
        )
    return cram_files
