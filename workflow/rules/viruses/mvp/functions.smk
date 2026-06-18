def get_sample_number(wildcards):
    """Return the 1-based row number of a sample_id in the MVP metadata.tsv"""
    return list(SAMPLES).index(wildcards.sample_id) + 1


def get_dram_input_fasta(wildcards):
    """Locate MVP_06's --DRAM fasta output (dynamic name, hence the checkpoint)"""
    checkpoints.viruses__mvp__06_functional_annotation.get()
    (names,) = glob_wildcards(
        VIR_MVP_FUNCTIONAL_ANNOTATION / "06_DRAM_V" / "{name}_DRAM_Input.fna"
    )
    return VIR_MVP_FUNCTIONAL_ANNOTATION / "06_DRAM_V" / f"{names[0]}_DRAM_Input.fna"


def get_dram_input_tsv(wildcards):
    """Locate MVP_06's --DRAM tsv output (dynamic name, hence the checkpoint)"""
    checkpoints.viruses__mvp__06_functional_annotation.get()
    (names,) = glob_wildcards(
        VIR_MVP_FUNCTIONAL_ANNOTATION / "06_DRAM_V" / "{name}_DRAM_Input.tsv"
    )
    return VIR_MVP_FUNCTIONAL_ANNOTATION / "06_DRAM_V" / f"{names[0]}_DRAM_Input.tsv"
