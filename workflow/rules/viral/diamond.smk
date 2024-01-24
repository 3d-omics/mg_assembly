rule _viral__diamond:
    input:
        proteins=GENOMAD / "{assembly_id}_summary" / "{assembly_id}_virus_proteins.faa",
        database=features["databases"]["diamond"],
    output:
        tsv=DIAMOND / "{assembly_id}.tsv",
    log:
        DIAMOND / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    shell:
        """
        diamond blastp \
            --threads {threads} \
            --db {input.database} \
            --query {input.proteins} \
            --outfmt 6 \
            --out {output.tsv} \
        2> {log} 1>&2
        """


rule viral__diamond:
    input:
        [DIAMOND / f"{assembly_id}.tsv" for assembly_id in ASSEMBLIES],
