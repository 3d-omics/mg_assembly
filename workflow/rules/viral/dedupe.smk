rule _viral__dedupe__unique_seqs:
    input:
        fastas=[CHECKV / f"{assembly_id}" / "all.fna" for assembly_id in ASSEMBLIES],
    output:
        fasta=DEDUPE / "dedupe.fa",
    log:
        DEDUPE / "dedupe.log",
    conda:
        "__environment__.yml"
    threads: 24
    params:
        fastas_comma=lambda w, input: ",".join(input),
        minimum_length=500,
    shell:
        """
        dedupe.sh \
            in={params.fastas_comma} \
            out={output.fasta} \
            minscaf={params.minimum_length} \
            overwrite=true \
            mergenames=t \
            exact=f \
            threads={threads} \
            usejni=t \
        2> {log} 1>&2
        """


rule viral__dedupe:
    input:
        DEDUPE / "dedupe.fa",
