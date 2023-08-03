def get_bams_for_concoct_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLY_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam",
        )
    return bam_files


def get_bais_for_concoct_binning(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = []
    for sample_id, library_id in samples_in_assembly:
        bam_files.append(
            ASSEMBLY_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam.bai",
        )
    return bam_files
