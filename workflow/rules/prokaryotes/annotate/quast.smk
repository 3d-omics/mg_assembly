use rule quast as prokaryotes__annotate__quast with:
    input:
        PROK_ANN / "drep.{secondary_ani}.fa.gz",
    output:
        directory(PROK_QUAST / "drep.{secondary_ani}"),
    log:
        PROK_QUAST / "drep.{secondary_ani}.log",


rule prokaryotes__annotate__quast__all:
    input:
        [PROK_QUAST / f"drep.{secondary_ani}" for secondary_ani in SECONDARY_ANIS],
