rule _viral__cluster__mmseqs__easy_cluster:
    input:
        fasta=DEDUPE / "dedupe.fa",
    output:
        fasta=MMSEQS / "cluster.fa",
    log:
        MMSEQS / "easy_cluster.log",
    conda:
        "__environment__.yml"
    params:
        prefix=MMSEQS / "results",
        tmpdir=MMSEQS,
        tmp_fasta=MMSEQS / "results_all_seqs.fasta",
    threads: 24
    shell:
        """
        mmseqs easy-cluster \
            {input.fasta} \
            {params.prefix} \
            {params.tmpdir} \
            --threads {threads} \
        2> {log} 1>&2

        mv {params.tmp_fasta} {output.fasta} 2>> {log} 1>&2
        """


rule viral__cluster__mmseqs:
    input:
        rules._viral__cluster__mmseqs__easy_cluster.output,
