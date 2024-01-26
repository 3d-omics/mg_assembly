rule _viral__mmseqs__easy_cluster:
    input:
        fasta=DEDUPE / "dedupe.fa",
    output:
        fasta=MMSEQS / "results_all_seqs.fasta",
    log:
        MMSEQS / "easy_cluster.log",
    conda:
        "__environment__.yml"
    params:
        prefix=MMSEQS / "results",
        tmpdir=MMSEQS,
    threads: 24
    shell:
        """
        mmseqs easy-cluster \
            {input.fasta} \
            {params.prefix} \
            {params.tmpdir} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule viral__mmseqs:
    input:
        rules._viral__mmseqs__easy_cluster.output,
