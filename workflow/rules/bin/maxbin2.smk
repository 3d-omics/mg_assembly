rule bin_maxbin2_prepare_one:
    """Compute coverages"""
    input:
        bams=get_bams_from_assembly_id,
    output:
        coverage=MAXBIN2 / "prepare" / "{assembly_id}.coverage",
    log:
        MAXBIN2 / "prepare" / "{assembly_id}.log",
    conda:
        "maxbin2.yml"
    shell:
        """
        (samtools coverage {input.bams} \
        | awk '{{print $1"\t"$5}}' \
        | grep -v '^#' \
        > {output.coverage} \
        ) 2> {log}
        """


rule bin_maxbin2_run_one:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=ASSEMBLE_RENAME / "{assembly_id}.fa",
        coverage=MAXBIN2 / "prepare" / "{assembly_id}.coverage",
    output:
        outdir=directory(MAXBIN2 / "bins" / "{assembly_id}/"),
    log:
        MAXBIN2 / "bins" / "{assembly_id}.log",
    conda:
        "maxbin2.yml"
    threads: 24
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
        """


rule bin_maxbin2:
    """Run MaxBin2 over all assemblies"""
    input:
        [MAXBIN2 / "bins" / f"{assembly_id}" for assembly_id in ASSEMBLIES],
