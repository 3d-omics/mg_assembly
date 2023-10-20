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
        rules.pre_eval_fastp_fastqc.input,
        rules.pre_eval_coverm.output,
        rules.pre_eval_samtools.input,
        rules.pre_eval_nonhost_fastqc.input,
        rules.pre_eval_kraken2.output,


rule pre_eval_with_singlem:
    input:
        rules.pre_eval.input,
        rules.pre_eval_singlem.output,


rule pre_eval_with_nonpareil:
    input:
        rules.pre_eval.input,
        rules.pre_eval_nonpareil.output,
