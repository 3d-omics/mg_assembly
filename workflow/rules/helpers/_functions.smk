def get_sample_and_library_from_assembly_id(assembly_id):
    """Get all the sample and library ids for a given assembly_id"""
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    return samples_in_assembly


def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024
