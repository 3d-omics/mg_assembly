rule preprocess__nonpareil__:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        cram=get_host_clean_cram,
    output:
        npa=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npa"),
        npc=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npc"),
        npl=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npl"),
        npo=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npo"),
    log:
        NONPAREIL / "run" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
        forward_fq=lambda w: NONPAREIL / "run" / f"{w.sample_id}.{w.library_id}_1.fq",
    shell:
        """
        samtools fastq \
            --threads {threads} \
            -1 {params.forward_fq} \
            -2 /dev/null \
            -0 /dev/null \
            -f 12 \
            {input.cram} \
        2> {log} 1>&2

        nonpareil \
            -s {params.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} \
        1>&2 || true

        rm --force --verbose {params.forward_fq} 2>> {log} 1>&2
        """


rule preprocess__nonpareil:
    """Aggregate all the nonpareil results into a single table"""
    input:
        [
            NONPAREIL / "run" / f"{sample_id}.{library_id}.{suffix}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for suffix in ["npa", "npc", "npl", "npo"]
        ],
    output:
        NONPAREIL / "nonpareil.tsv",
    log:
        NONPAREIL / "nonpareil.log",
    conda:
        "__environment__.yml"
    params:
        input_dir=NONPAREIL / "run",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """
