# Output naming (PHROGS/PFAM/DRAM-v inputs) is dynamic, so this is a checkpoint;
# see get_dram_input_fasta/get_dram_input_tsv in mvp/functions.smk.


checkpoint viruses__mvp__06_functional_annotation:
    """Functional annotation (PHROGS/PFAM) of the representative catalog + DRAM-v input prep"""
    input:
        metadata=VIR_MVP / "metadata.tsv",
        representative_fasta=VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences.fna",
    output:
        directory(VIR_MVP_FUNCTIONAL_ANNOTATION),
    log:
        VIR_MVP / "mvp_06_functional_annotation.log",
    conda:
        ENVS / "mvp.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=24 * 60,
    params:
        fasta_files=params["viral"]["mvp"]["functional_annotation"]["fasta_files"],
        phrogs_evalue=params["viral"]["mvp"]["functional_annotation"]["phrogs_evalue"],
        phrogs_score=params["viral"]["mvp"]["functional_annotation"]["phrogs_score"],
        pfam_evalue=params["viral"]["mvp"]["functional_annotation"]["pfam_evalue"],
        pfam_score=params["viral"]["mvp"]["functional_annotation"]["pfam_score"],
    shell:
        """
        mvip MVP_06_do_functional_annotation \
            --input {VIR_MVP} \
            --metadata {input.metadata} \
            --fasta_files {params.fasta_files} \
            --PHROGS_evalue {params.phrogs_evalue} \
            --PHROGS_score {params.phrogs_score} \
            --PFAM_evalue {params.pfam_evalue} \
            --PFAM_score {params.pfam_score} \
            --DRAM \
            --threads {threads} \
        2> {log} 1>&2
        """
