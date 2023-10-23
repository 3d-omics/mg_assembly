include: "functions.smk"
include: "fastp.smk"
include: "bowtie2.smk"
include: "fastqc.smk"
include: "nonpareil.smk"
include: "singlem.smk"
include: "coverm.smk"
include: "kraken2.smk"
include: "samtools.smk"


# FASTP report files ----
rule pre_eval_fastp:
    input:
        [
            FASTP / f"{sample_id}.{library_id}_fastp.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],


rule pre_eval:
    input:
        rules.pre_fastp_fastqc.input,
        rules.pre_coverm.output,
        rules.pre_samtools.input,
        rules.pre_nonhost_fastqc.input,
        rules.pre_kraken2.output,


rule pre_eval_with_singlem:
    input:
        rules.pre_eval.input,
        rules.pre_singlem.output,


rule pre_eval_with_nonpareil:
    input:
        rules.pre_eval.input,
        rules.pre_nonpareil.output,


rule pre_run:
    input:
        rules.pre_bowtie2_extract_nonhost_all.input,


rule pre:
    input:
        rules.pre_run.input,
        rules.pre_eval.input,
