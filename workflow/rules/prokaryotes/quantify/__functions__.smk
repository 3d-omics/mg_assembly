def _compose_tsv_for_coverm(wildcards, mode):
    """Auxiliary function to compose coverm files"""
    assert mode in ["contig", "genome"]
    method = wildcards.method
    secondary_ani = wildcards.secondary_ani
    return [
        COVERM / mode / method / f"drep.{secondary_ani}" / f"{sample_id}.{library_id}.tsv.gz"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]


def get_tsvs_for_dereplicate_coverm_genome(wildcards):
    """Get the tsv files for the coverm genome rule."""
    return _compose_tsv_for_coverm(wildcards, mode="genome")


def get_tsvs_for_dereplicate_coverm_contig(wildcards):
    """Get the tsv files for the coverm contig rule."""
    return _compose_tsv_for_coverm(wildcards, mode="contig")
