def get_forwards_for_megahit(wildcards):
    """Get the forward files for megahit"""
    sample_id = wildcards.sample_id
    library_ids = samples[samples.sample_id == sample_id].library_id.tolist()
    forward_filenames = [
        NONHOST / f"{sample_id}.{library_id}_1.fq.gz" for library_id in library_ids
    ]
    return forward_filenames


def get_reverses_for_megahit(wildcards):
    """Get the reverse files for megahit"""
    sample_id = wildcards.sample_id
    library_ids = samples[samples.sample_id == sample_id].library_id.tolist()
    reverse_filenames = [
        NONHOST / f"{sample_id}.{library_id}_2.fq.gz" for library_id in library_ids
    ]
    return reverse_filenames


def aggregate_forwards_for_megahit(wildcards):
    """Put all the forwards together separated by a comma"""
    return ",".join([str(x) for x in get_forwards_for_megahit(wildcards)])


def aggregate_reverses_for_megahit(wildcards):
    """Put all the reverses together separated by a comma"""
    return ",".join([str(x) for x in get_reverses_for_megahit(wildcards)])
