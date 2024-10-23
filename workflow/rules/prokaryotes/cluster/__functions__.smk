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
    """Given an assembly_id, get all the bam files for that assembly."""
    return [f"{bam}.bai" for bam in get_bams_from_assembly_id(wildcards)]


def compose_bams_for_metabat2_run(wildcards):
    """Given an assemblu_id, get all the bam files that will be generated in metabat2"""
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = [
        METABAT2 / f"{assembly_id}.{sample_id}.{library_id}.bam"
        for sample_id, library_id in samples_in_assembly
    ]
    return bam_files
