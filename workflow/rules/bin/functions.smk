def get_bams_from_assembly_id(wildcards):
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    bam_files = [
        ASSEMBLE_BOWTIE2 / f"{assembly_id}.{sample_id}.{library_id}.bam"
        for sample_id, library_id in samples_in_assembly
    ]
    return bam_files


def get_bais_from_assembly_id(wildcards):
    bams = get_bams_from_assembly_id(wildcards)
    return [f"{bam}.bai" for bam in bams]


# Magscot ----
def compose_out_prefix_for_bin_magscot_run_one(wildcards):
    return MAGSCOT / wildcards.assembly_id / "magscot"


# Concoct ----
def compose_basename_for_concoct_run_one(wildcards):
    return CONCOCT / "run" / wildcards.assembly_id


# MaxBin2 ----
def compose_out_prefix_for_maxbin2_run_one(wildcards):
    return MAXBIN2 / "bins" / wildcards.assembly_id


# MetaBat2 ---
def compose_bins_prefix_for_metabat2_run_one(wildcards):
    return METABAT2 / "bins" / wildcards.assembly_id / wildcards.assembly_id
