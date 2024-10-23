def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024


def get_attempt(wildcards, attempt):
    """Get the number of attempt in resources"""
    return attempt


def compose_rg_id(w):
    """Compose the read group ID for bowtie2"""
    return f"{w.sample_id}_{w.library_id}"


def compose_rg_extra(w):
    """Compose the read group extra information for bowtie2"""
    return f"LB:truseq_{w.library_id}\\tPL:Illumina\\tSM:{w.sample_id}"
