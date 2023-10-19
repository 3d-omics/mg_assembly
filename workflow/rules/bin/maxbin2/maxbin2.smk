include: "functions.smk"
include: "run.smk"


rule maxbin2:
    input:
        [MAXBIN2 / "bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
