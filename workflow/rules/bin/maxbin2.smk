rule _bin__maxbin2__prepare:
    """Compute coverages"""
    input:
        bams=get_bams_from_assembly_id,
    output:
        coverage=MAXBIN2 / "prepare" / "{assembly_id}.coverage",
    log:
        MAXBIN2 / "prepare" / "{assembly_id}.log",
    conda:
        "maxbin2.yml"
    resources:
        runtime=6 * 60,
        mem_mb=8 * 1024,
    shell:
        """
        ( samtools coverage {input.bams} \
        | awk '{{print $1"\t"$5}}' \
        | grep -v '^#' \
        > {output.coverage} \
        ) 2> {log}
        """


rule _bin__maxbin2__run:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
        coverage=MAXBIN2 / "prepare" / "{assembly_id}.coverage",
    output:
        outdir=directory(MAXBIN2 / "bins" / "{assembly_id}"),
    log:
        MAXBIN2 / "bins" / "{assembly_id}.log",
    conda:
        "maxbin2.yml"
    threads: 4
    params:
        seed=1,
        out_prefix=compose_out_prefix_for_maxbin2_run_one,
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
    shell:
        """
        mkdir --parents {output.outdir}

        run_MaxBin.pl \
            -thread {threads} \
            -contig {input.assembly} \
            -out {output.outdir}/maxbin2 \
            -abund {input.coverage} \
        2> {log} 1>&2

        rename 's/\.fasta$/.fa/' {output.outdir}/*.fasta 2>> {log}
        """


rule bin__maxbin2:
    """Run MaxBin2 over all assemblies"""
    input:
        [MAXBIN2 / "bins" / assembly_id for assembly_id in ASSEMBLIES],
