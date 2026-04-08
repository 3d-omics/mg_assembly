rule viruses__cluster__mmseqs__easy_cluster:
    input:
        fasta=VIR_CLUSTER / "bbmap.clean.fa.gz",
    output:
        all_seq=VIR_CLUSTER / "mmseqs.all_seqs.fa.gz",
        cluster=VIR_CLUSTER / "mmseqs.cluster.tsv.gz",
        rep_seq=VIR_CLUSTER / "mmseqs.rep_seq.fa.gz",
    log:
        VIR_CLUSTER / "mmseqs.easy_cluster.log",
    conda:
        ENVS / "mmseqs.yml"
    params:
        prefix=VIR_CLUSTER / "tmp",
        tmpdir=VIR_CLUSTER,
    shadow:
        "minimal"
    threads: 24
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


rule viruses__cluster__mmseqs__all:
    input:
        rules.viruses__cluster__mmseqs__easy_cluster.output,
