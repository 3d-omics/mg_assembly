def get_sample_and_library_from_assembly_id(assembly_id):
    """Get all the sample and library ids for a given assembly_id"""
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    return samples_in_assembly


def _get_reads_from_assembly_id(wildcards, end):
    """Get the file_end for megahit"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    return [
        PRE_BOWTIE2 / f"{sample_id}.{library_id}_{end}.fq.gz"
        for sample_id, library_id in samples_in_assembly
    ]


def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return _get_reads_from_assembly_id(wildcards, end="forward")


def get_reverses_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return _get_reads_from_assembly_id(wildcards, end="reverse")


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = [
        str(forward_)
        for forward_ in _get_reads_from_assembly_id(wildcards, end="forward")
    ]
    return ",".join(forwards)


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = [
        str(reverse_)
        for reverse_ in _get_reads_from_assembly_id(wildcards, end="reverse")
    ]
    return ",".join(reverses)
