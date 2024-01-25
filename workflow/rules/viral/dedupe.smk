rule _viral__dedupe__unique_seqs:
    input:
        fastas = [CHECKV / f"{assembly_id}" / "all.fna" for assembly_id in ASSEMBLIES],
    output:
        fasta=DEDUPE / "unique_seqs.fa",
        stats=DEDUPE / "stats.tsv"
    log:
        DEDUPE / "unique_seqs.log",
    conda:
        "__environment__.yml"
    threads:
        24
    params:
        fastas=",".join([str(CHECKV / f"{assembly_id}" / "all.fna") for assembly_id in ASSEMBLIES]),
        workdir=DEDUPE,
        min_scaf=500,
    shell:
        """
        dedupe.sh \
            in={params.fastas} \
            out={params.workdir} \
            csf={output.stats} \
            minscaf={params.min_scaf} \
            mergenames=t \
            exact=f \
            threads={threads} \
            usejni=t \
        2> {log} 1>&2
        """


rule viral__dedupe:
    input:
        DEDUPE / "unique_seqs.fa"
