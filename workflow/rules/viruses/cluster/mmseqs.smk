rule viruses__cluster__mmseqs:
    input:
        fasta=VIR_DEDUPE / "clean.fa.gz",
    output:
        all_seq=VIR_MMSEQS / "all_seqs.fa.gz",
        cluster=VIR_MMSEQS / "cluster.tsv.gz",
        rep_seq=VIR_MMSEQS / "rep_seq.fa.gz",
    log:
        VIR_MMSEQS / "easy_cluster.log",
    conda:
        "../../../environments/mmseqs.yml"
    params:
        prefix=VIR_MMSEQS / "tmp",
        tmpdir=VIR_MMSEQS,
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
        rules.viruses__cluster__mmseqs.output,
