def get_bams_for_metawrap_metawrap_prepare(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files


def compose_metawrap_working_folder(wildcards):
    assembly_id = wildcards.assembly_id
    completeness = params["metabin"]["metawrap"]["refinement"]["completeness"]
    contamination = params["metabin"]["metawrap"]["refinement"]["contamination"]
    folder = (
        METAWRAP_REFINEMENT
        / f"{assembly_id}/metawrap_{completeness}_{contamination}_bins"
    )
    return folder
