rule viruses__quantify__index:
    """Index dereplicader"""
    input:
        contigs=MMSEQS / "rep_seq.fasta",
    output:
        mock=touch(VINDEX / "viruses"),
    log:
        VINDEX / "virues.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["quantify"]["bowtie2-build"]["memory_gb"]),
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.contigs} \
            {output.mock} \
        2> {log} 1>&2
        """
