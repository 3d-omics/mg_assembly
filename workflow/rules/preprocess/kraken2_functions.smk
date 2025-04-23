def get_libraries_from_sample(wildcards):
    """
    Get the sample_id, library_id from the associated sample_id
    """

    sample_library_filt = [
        [sample_id, library_id] for sample_id, library_id in SAMPLE_LIBRARY 
        if sample_id == wildcards.sample_id
    ]

    return sample_library_filt