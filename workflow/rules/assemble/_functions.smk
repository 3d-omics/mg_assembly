# Megahit
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


def compose_out_dir_for_assemble_megahit_one(wildcards):
    """Compose output folder"""
    return MEGAHIT / wildcards.assembly_id


# Concoct, maxbin and metabat need bams to get coverage
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


def compose_basename_for_concoct_run_one(wildcards):
    """Compose the basename for the concoct run"""
    return CONCOCT / "run" / wildcards.assembly_id


def compose_out_prefix_for_maxbin2_run_one(wildcards):
    """Compose the output folder for maxbin2"""
    return MAXBIN2 / "bins" / wildcards.assembly_id


def compose_bins_prefix_for_metabat2_run_one(wildcards):
    """Compose the output folder for metabat2"""
    return METABAT2 / "bins" / wildcards.assembly_id / wildcards.assembly_id


# Magscot ----
def compose_out_prefix_for_bin_magscot_run_one(wildcards):
    """Compose the output folder for magscot"""
    return MAGSCOT / wildcards.assembly_id / "magscot"
