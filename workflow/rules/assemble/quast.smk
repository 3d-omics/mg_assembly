use rule quast as assemble__quast with:
    input:
        ASMB_MEGAHIT / "{assembly_id}.fa.gz",
    output:
        directory(ASMB_QUAST / "{assembly_id}"),
    log:
        ASMB_QUAST / "{assembly_id}.log",


rule assemble__quast__all:
    input:
        [ASMB_QUAST / f"{assembly_id}" for assembly_id in ASSEMBLIES],
