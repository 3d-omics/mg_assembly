rule reads_eval_fastqc_one:
    input:
        fq=READS / "{sample}.{library}_{end}.fq.gz",
    output:
        html=READS / "{sample}.{library}_{end}_fastqc.html",
        zip_=READS / "{sample}.{library}_{end}_fastqc.zip",
    log:
        READS / "{sample}.{library}_{end}_fastqc.log",
    conda:
        "../envs/reads.yml"
    shell:
        """
        fastqc --outdir {READS} {input.fq} &> {log}
        """


rule reads_eval:
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],
