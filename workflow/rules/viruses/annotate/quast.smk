use rule quast as viruses__annotate__quast with:
    input:
        VIR_CLUSTER / "mmseqs.rep_seq.fa.gz",
    output:
        directory(VIR_QUAST),
    log:
        VIR / "quast.log",


rule viruses__annotate__quast__all:
    input:
        rules.viruses__annotate__quast.output,
