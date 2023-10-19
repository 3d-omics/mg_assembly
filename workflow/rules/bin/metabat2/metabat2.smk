include: "functions.smk"
include: "run.smk"


rule metabat2:
    input:
        [METABAT2 / f"bins/{assembly_id}" for assembly_id in ASSEMBLIES],
