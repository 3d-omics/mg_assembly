rule _viral__mmseqs__easy_cluster:
    input:
        fasta=DEDUPE / "dedupe.fa",
    output:
        MMSEQS / "easy_cluster.fa"
    log:
        MMSEQS / "easy_cluster.log"
    conda:
        "__environment__.yml"
    params:
        prefix=
        tmpdir=
    threads:
        24
    shell:
        """
        mmseqs easy-cluster \
            {input.fasta} \
            {params.prefix} \
            {params.tmpdir} \
            --threads {threads} \
        2> {log} 1>&2
        """
