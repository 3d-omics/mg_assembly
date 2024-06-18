rule preprocess__fastp__:
    """Run fastp on one PE library"""
    input:
        forward_=READS / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=READS / "{sample_id}.{library_id}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample_id}.{library_id}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample_id}.{library_id}_2.fq.gz"),
        html=FASTP / "{sample_id}.{library_id}_fastp.html",
        json=FASTP / "{sample_id}.{library_id}_fastp.json",
    log:
        FASTP / "{sample_id}.{library_id}.log",
    params:
        length_required=params["preprocess"]["fastp"]["length_required"],
        forward_adapter=get_forward_adapter,
        reverse_adapter=get_reverse_adapter,
    conda:
        "__environment__.yml"
    shell:
        """
        fastp \
            --in1 <(gzip --decompress --stdout {input.forward_}) \
            --in2 <(gzip --decompress --stdout {input.reverse_}) \
            --out1 >(pigz --processes {threads} --fast > {output.forward_}) \
            --out2 >(pigz --processes {threads} --fast > {output.reverse_}) \
            --adapter_sequence {params.forward_adapter} \
            --adapter_sequence_r2 {params.reverse_adapter} \
            --html {output.html} \
            --json {output.json} \
            --length_required {params.length_required} \
            --thread {threads} \
            --trim_poly_g \
            --trim_poly_x \
            --verbose \
        2> {log} 1>&2
        """


rule preprocess__fastp__import__:
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
    output:
        cram=FASTP / "{sample_id}.{library_id}.cram",
    log:
        FASTP / "{sample_id}.{library_id}.cram.log",
    conda:
        "__environment__.yml"
    group:
        "sample"
    shell:
        """
        samtools import \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --output-fmt CRAM \
            --threads {threads} \
            -o {output.cram} \
        2> {log} 1>&2
        """


rule preprocess__fastp:
    """Run fastp over all libraries"""
    input:
        cram=[
            FASTP / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        json=[
            FASTP / f"{sample_id}.{library_id}_fastp.json"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
