def get_tsvs_for_dereplicate_coverm_genome(wildcards):
    """Get the tsv files for the coverm genome rule."""
    method = wildcards.method
    tsv_files = [
        COVERM / "genome" / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def get_tsvs_for_dereplicate_coverm_contig(wildcards):
    """Get the tsv files for the coverm contig rule."""
    method = wildcards.method
    tsv_files = [
        COVERM / "contig" / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
    return tsv_files


def compose_input_dir_for_dereplicate_coverm_contig_method(wildcards):
    """Get the input folder for coverm contig"""
    return COVERM / "contig" / wildcards.method


def compose_input_dir_for_dereplicate_coverm_genome_method(wildcards):
    """Get the input folder for coverm genome"""
    return COVERM / "genome" / wildcards.method
