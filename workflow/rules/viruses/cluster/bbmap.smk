rule viruses__cluster__bbmap__dedupe:
    input:
        VIR_CLUSTER / "genomad_virus.fna.gz"
    output:
        fasta=VIR_CLUSTER / "bbmap.dedupe.fa.gz",
    log:
        VIR_CLUSTER / "bbmap.dedupe.log",
    conda:
        "../../../environments/bbmap.yml"
    params:
        fastas_comma=lambda w, input: ",".join(input),
        minimum_length=500,
    threads: 24
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


rule viruses__cluster__bbmap__clean:
    """
    Clean up the deduped fasta file since merged sequences headers contain multiple ">"
    """
    input:
        VIR_CLUSTER / "bbmap.dedupe.fa.gz",
    output:
        VIR_CLUSTER / "bbmap.clean.fa.gz",
    log:
        VIR_CLUSTER / "bbmap.clean.log",
    conda:
        "../../../environments/bbmap.yml"
    shell:
        """
        ( seqtk seq {input} \
        | cut --fields 1,2 -d ">" \
        | tr "|" "_" \
        | bgzip \
        > {output} \
        ) 2> {log}
        """


rule viruses__cluster__bbmap__all:
    input:
        rules.viruses__cluster__bbmap__clean.output,
