# fastp ----
def get_adapter(wildcards, forward_or_reverse):
    """Get forward or reverse adapter"""
    assert forward_or_reverse in ["forward_adapter", "reverse_adapter"]
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ][forward_or_reverse].tolist()[0]


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return get_adapter(wildcards, "forward_adapter")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return get_adapter(wildcards, "reverse_adapter")


# bowtie2 ----
def compose_rg_id(wildcards):
    """Compose the read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose the read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.sample_id}"


def get_input_cram_for_host_mapping(wildcards):
    """Get the input cram file for host mapping"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    genome = wildcards.genome
    if genome == HOST_NAMES[0]:
        return FASTP / f"{sample_id}.{library_id}.cram"
    genome_index = HOST_NAMES.index(genome)
    previous_genome = HOST_NAMES[genome_index - 1]
    return PRE_BOWTIE2 / f"{previous_genome}" / f"{sample_id}.{library_id}.cram"


# clean ----
def get_host_clean_cram(wildcards):
    """Get the input file that is clean from hosts"""
    last_genome = HOST_NAMES[-1]
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample_id}.{library_id}.cram"
    return PRE_BOWTIE2 / last_genome / f"{sample_id}.{library_id}.cram"
