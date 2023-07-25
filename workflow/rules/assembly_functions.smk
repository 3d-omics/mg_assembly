def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    assembly_id = wildcards.assembly_id
    forward_filenames = samples[
        samples.assembly_id == assembly_id
    ].forward_filename.values.tolist()
    return forward_filenames


def get_reverses_from_assembly_id(wildcards):
    """Get the reverse files for megahit"""
    assembly_id = wildcards.assembly_id
    reverse_filenames = samples[
        samples.assembly_id == assembly_id
    ].reverse_filename.values.tolist()
    return reverse_filenames


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = get_forwards_from_assembly_id(wildcards)
    return ",".join(forwards)


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = get_reverses_from_assembly_id(wildcards)
    return ",".join(reverses)
