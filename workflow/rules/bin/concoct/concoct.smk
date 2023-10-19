include: "functions.smk"
include: "run.smk"


rule concoct:
    input:
        [CONCOCT / "fasta_bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
