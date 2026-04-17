rule prokaryotes__cluster__maxbin2:
    """Run MaxBin2 over a single assembly"""
    input:
        assembly=ASMB_MEGAHIT / "{assembly_id}.fa.gz",
        bams=get_bams_from_assembly_id,
    output:
        workdir=directory(PROK_MAXBIN2 / "{assembly_id}"),
    log:
        PROK_MAXBIN2 / "{assembly_id}.log",
    conda:
        ENVS / "maxbin2.yml"
    threads: 4
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=24 * 60,
    params:
        coverage=lambda w: PROK_MAXBIN2 / w.assembly_id / "maxbin2.coverage",
    shell:
        """
        mkdir --parents {output.workdir}

        ( samtools coverage {input.bams} \
        | awk '{{print $1"\\t"$5}}' \
        | grep -v '^#' \
        ) > {params.coverage} \
        2> {log}

        run_MaxBin.pl \
            -thread {threads} \
            -contig {input.assembly} \
            -out {output.workdir}/maxbin2 \
            -abund {params.coverage} \
        2> {log} 1>&2

        rename \
            's/\\.fasta$/.fa/' \
            {output.workdir}/*.fasta \
        2>> {log}

        bgzip \
            {output.workdir}/*.fa \
        2>> {log} 1>&2

        rm \
            --recursive \
            --force \
            {output.workdir}/maxbin.{{coverage,log,marker,noclass,summary,tooshort}} \
            {output.workdir}/maxbin2.marker_of_each_bin.tar.gz \
        2>> {log} 1>&2
        """


rule prokaryotes__cluster__maxbin2__all:
    """Run MaxBin2 over all assemblies"""
    input:
        [PROK_MAXBIN2 / assembly_id for assembly_id in ASSEMBLIES],
