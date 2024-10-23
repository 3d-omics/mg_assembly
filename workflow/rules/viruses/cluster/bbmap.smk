rule viruses__cluster__bbmap__dedupe:
    input:
        fastas=[GENOMADC / f"{assembly_id}_virus.fna.gz" for assembly_id in ASSEMBLIES],
    output:
        fasta=DEDUPE / "dedupe.fa.gz",
    log:
        DEDUPE / "bbmap.log",
    conda:
        "../../../environments/bbmap.yml"
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


rule viruses__cluster__bbmap__all:
    input:
        rules.viruses__cluster__bbmap__dedupe.output,
