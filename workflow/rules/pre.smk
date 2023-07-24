rule pre_fastp_trim_one:
    """Run fastp on one library"""
    input:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample}.{library}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample}.{library}_2.fq.gz"),
        unpaired1=temp(FASTP / "{sample}.{library}_u1.fq.gz"),
        unpaired2=temp(FASTP / "{sample}.{library}_u2.fq.gz"),
        html=FASTP / "{sample}.{library}_fastp.html",
        json=FASTP / "{sample}.{library}_fastp.json",
    log:
        FASTP / "{sample}.{library}.log",
    conda:
        "../envs/pre.yml"
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["fastp"]["extra"],
        length_required=params["fastp"]["length_required"],
    threads: 16
    resources:
        mem_mb=4 * 1024,
        runtime=240,
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 >(gzip --fast > {output.forward_}) \
            --out2 >(gzip --fast > {output.reverse_}) \
            --unpaired1 >(gzip --fast > {output.unpaired1}) \
            --unpaired2 >(gzip --fast > {output.unpaired2}) \
            --html {output.html} \
            --json {output.json} \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule pre_fastp_trim_all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule pre_fastp_fastqc_one:
    """Run fastqc on one library from fastp output"""
    input:
        fq=FASTP / "{sample}.{library}_{end}.fq.gz",
    output:
        html=FASTP / "{sample}.{library}_{end}_fastqc.html",
        zip_=FASTP / "{sample}.{library}_{end}_fastqc.zip",
    log:
        FASTP / "{sample}.{library}_{end}_fastqc.log",
    conda:
        "../envs/pre.yml"
    shell:
        """
        fastqc \
            --outdir {FASTP} \
            --threads 1 \
            {input.fq} \
        2> {log} 1>&2
        """


rule pre_fastp_fastqc_all:
    """Run fastqc over all libraries after fastp"""
    input:
        [
            FASTP / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
            for extension in "html zip".split(" ")
        ],


rule pre:
    input:
        rules.pre_fastp_trim_all.input,
        rules.pre_fastp_fastqc_all.input,
