def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
    ]["reverse_adapter"].tolist()[0]


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
    return f"{lb_field}\t" + f"PL:Illumina\t" + f"{sm_field}"


def compose_prefix_for_nonpareil(wildcards):
    """Compose prefix for nonpareil output files"""
    return NONPAREIL / f"{wildcards.sample_id}.{wildcards.library_id}"
