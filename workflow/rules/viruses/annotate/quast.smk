use rule quast as viruses__annotate__quast with:
    input:
        VIR_MVP_CLUSTERING
        / "MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences.fna",
    output:
        directory(VIR_QUAST),
    log:
        VIR / "quast.log",


rule viruses__annotate__quast__all:
    input:
        rules.viruses__annotate__quast.output,
