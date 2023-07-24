def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    return samples[(samples["sample_id"] == wildcards.sample)][
        "forward_filename"
    ].tolist()[0]


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    return samples[(samples["sample_id"] == wildcards.sample)][
        "reverse_filename"
    ].tolist()[0]
