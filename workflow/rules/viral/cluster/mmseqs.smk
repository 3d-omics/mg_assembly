rule _viral__cluster__mmseqs__easy_cluster:
    input:
        fasta=DEDUPE / "dedupe.fa",
    output:
        all_seq=MMSEQS / "all_seqs.fasta",
        cluster=MMSEQS / "cluster.tsv",
        rep_seq=MMSEQS / "rep_seq.fasta",
    log:
        MMSEQS / "easy_cluster.log",
    conda:
        "__environment__.yml"
    params:
        prefix=MMSEQS / "tmp",
        tmpdir=MMSEQS,
    threads: 24
    shadow: "minimal",
    shell:
        """
        mmseqs easy-cluster \
            {input.fasta} \
            {params.prefix} \
            {params.tmpdir} \
            --threads {threads} \
        2> {log} 1>&2

        mv {params.tmpdir}/tmp_all_seqs.fasta {output.all_seq} 2>> {log} 1>&2
        mv {params.tmpdir}/tmp_cluster.tsv {output.cluster} 2>> {log} 1>&2
        mv {params.tmpdir}/tmp_rep_seq.fasta {output.rep_seq} 2>> {log} 1>&2
        """


rule viral__cluster__mmseqs:
    input:
        rules._viral__cluster__mmseqs__easy_cluster.output,
