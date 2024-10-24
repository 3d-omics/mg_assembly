rule viruses__cluster__mmseqs__:
    input:
        fasta=DEDUPE / "dedupe.fa.gz",
    output:
        all_seq=MMSEQS / "all_seqs.fasta.gz",
        cluster=MMSEQS / "cluster.tsv.gz",
        rep_seq=MMSEQS / "rep_seq.fasta.gz",
    log:
        MMSEQS / "easy_cluster.log",
    conda:
        "__environment__.yml"
    params:
        prefix=MMSEQS / "tmp",
        tmpdir=MMSEQS,
    shadow:
        "minimal"
    shell:
        """
        mmseqs easy-cluster \
            {input.fasta} \
            {params.prefix} \
            {params.tmpdir} \
            --threads {threads} \
        2> {log} 1>&2

        bgzip \
            --threads {threads} \
            --stdout \
            {params.tmpdir}/tmp_all_seqs.fasta \
        > {output.all_seq} \
        2>> {log}

        bgzip \
            --threads {threads} \
            --stdout \
            {params.tmpdir}/tmp_cluster.tsv \
        > {output.cluster} \
        2>> {log}

        bgzip \
            --threads {threads} \
            --stdout \
            {params.tmpdir}/tmp_rep_seq.fasta \
        > {output.rep_seq} \
        2>> {log}
        """


rule viruses__cluster__mmseqs:
    input:
        rules.viruses__cluster__mmseqs__.output,
