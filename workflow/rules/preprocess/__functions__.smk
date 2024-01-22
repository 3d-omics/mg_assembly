# fastp
def get_adapter(wildcards, end):
    """Get the adapter of the en from a file"""
    assert end in ["forward", "reverse"]
    end = "forward_adapter" if end == "forward" else "reverse_adapter"
    adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][end].tolist()[0]
    return adapter

def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return get_adapter(wildcards, end="forward")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return get_adapter(wildcards, end="reverse")


# bowtie2
def compose_rg_id(wildcards):
    """Compose read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose read group extra information for bowtie2"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    lb_field = f"LB:{sample_id}_{library_id}"
    pl_field = "PL:Illumina"
    sm_field = f"SM:{sample_id}"
    return f"{lb_field}\t{pl_field}\t{sm_field}"


def get_input_file_for_host_mapping(wildcards, end):
    """Compose the input file for host mapping"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample_id}.{wildcards.library_id}_{end}.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        PRE_BOWTIE2
        / f"non{prev_genome}"
        / f"{wildcards.sample_id}.{wildcards.library_id}_{end}.fq.gz"
    )


def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
   return get_input_file_for_host_mapping(wildcards, end = "forward")


def get_input_reverse_for_host_mapping(wildcards):
    """Compose the forward input file"""
   return get_input_file_for_host_mapping(wildcards, end = "reverse")


# nonpareil
def compose_prefix_for_nonpareil(wildcards):
    """Compose prefix for nonpareil output files"""
    return NONPAREIL / f"{wildcards.sample_id}.{wildcards.library_id}"


def get_cram_for_preprocess_eval_cram_to_mapped_bam(wildcards):
    """Get the cram file of the last host"""
    genome = LAST_HOST
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return PRE_BOWTIE2 / genome / f"{sample_id}.{library_id}.cram"


def get_crai_for_preprocess_eval_cram_to_mapped_bam(wildcards):
    """Get the crai file of the last host"""
    bam = get_cram_for_pre_eval_cram_to_mapped_bam(wildcards)
    return f"{bam}.bai"

# kraken2
def get_kraken2_database(wildcards):
    """Get kraken2 database path from the name"""
    return features["databases"]["kraken2"][wildcards.kraken_db]


def compose_out_folder_for_eval_kraken2_assign_all(wildcards):
    """Just compose the output folder"""
    return KRAKEN2 / wildcards.kraken_db


# finals
def get_final_forward_from_pre(wildcards):
    """Get the last host forward file or the result from FASTP"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    last_host = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample_id}.{library_id}_1.fq.gz"
    return PRE_BOWTIE2 / f"non{last_host}" / f"{sample_id}.{library_id}_1.fq.gz"


def get_final_reverse_from_pre(wildcards):
    """Get the last host reverse file or the result from FASTP"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    last_host = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample_id}.{library_id}_2.fq.gz"
    return PRE_BOWTIE2 / f"non{last_host}" / f"{sample_id}.{library_id}_2.fq.gz"
