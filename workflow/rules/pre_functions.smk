def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]
