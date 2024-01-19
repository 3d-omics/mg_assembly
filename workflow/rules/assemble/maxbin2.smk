rule _assemble__maxbin2__run:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=MEGAHIT / "{assembly_id}.fa.gz",
        crams=get_crams_from_assembly_id,
    output:
        outdir=directory(MAXBIN2 / "{assembly_id}"),
    log:
        MAXBIN2 / "{assembly_id}.log",
    conda:
        "__environment__.yml"
    threads: 4
    params:
        seed=1,
        out_prefix=lambda w: MAXBIN2 / f"{w.assembly_id}",
        coverage=lambda w: MAXBIN2 / f"{w.assembly_id}/maxbin2.coverage",
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    shell:
        """
        mkdir --parents {output.outdir}

        ( samtools coverage {input.crams} \
        | awk '{{print $1"\\t"$5}}' \
        | grep -v '^#' \
        ) > {params.coverage} \
        2> {log}

        run_MaxBin.pl \
            -thread {threads} \
            -contig {input.assembly} \
            -out {output.outdir}/maxbin2 \
            -abund {params.coverage} \
        2> {log} 1>&2

        rename \
            's/\\.fasta$/.fa/' \
            {output.outdir}/*.fasta \
        2>> {log}

        find \
            {output.outdir} \
            -name "*.fa" \
            -exec pigz --best --verbose {{}} \; \
        2>> {log} 1>&2
        """


rule assemble__maxbin2:
    """Run MaxBin2 over all assemblies"""
    input:
        [MAXBIN2 / assembly_id for assembly_id in ASSEMBLIES],
