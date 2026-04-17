include: "kraken2_functions.smk"


rule preprocess__kraken2__join_libraries:
    """Join all libraries for a single sample"""
    input:
        forwards=lambda w: [
            PRE_FASTP / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in get_libraries_from_sample(w)
        ],
        reverses=lambda w: [
            PRE_FASTP / f"{sample_id}.{library_id}_2.fq.gz"
            for sample_id, library_id in get_libraries_from_sample(w)
        ],
    output:
        forwards=temp(PRE_KRAKEN2 / "samples" / "{sample_id}_1.fq.gz"),
        reverses=temp(PRE_KRAKEN2 / "samples" / "{sample_id}_2.fq.gz"),
    log:
        PRE_KRAKEN2 / "samples" / "{sample_id}.log",
    conda:
        ENVS / "kraken2.yml"
    shell:
        """
        cat {input.forwards} > {output.forwards} 2> {log}
        cat {input.reverses} > {output.reverses} 2>> {log}
        """


rule preprocess__kraken2__assign:
    """
Run kraken2 over all samples at once using the /dev/shm/ trick.

NOTE:
    - /dev/shm may be not empty after the job is done.
    - Specify twice the amount of RAM needed: Linux systems usually
        come configured with /dev/shm to be half the RAM size
    - After read classification the report generation step uses suddenly
        ~10GB of RAM per sample processed in parallel. The 2x RAM usage
        comes handy to avoid OOM errors.
"""
    input:
        forwards=[
            PRE_KRAKEN2 / "samples" / f"{sample_id}_1.fq.gz" for sample_id in SAMPLES
        ],
        rerverses=[
            PRE_KRAKEN2 / "samples" / f"{sample_id}_2.fq.gz" for sample_id in SAMPLES
        ],
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
    output:
        out_gzs=[
            PRE_KRAKEN2 / "{kraken2_db}" / f"{sample_id}.out.gz"
            for sample_id in SAMPLES
        ],
        reports=[
            PRE_KRAKEN2 / "{kraken2_db}" / f"{sample_id}.k2report"
            for sample_id in SAMPLES
        ],
    log:
        PRE_KRAKEN2 / "{kraken2_db}.log",
    conda:
        ENVS / "kraken2.yml"
    threads: 24
    resources:
        mem_mb=2 * 800 * 1024,  # Use twice the size of the database, we use /dev/shm
        runtime=24 * 60,
    params:
        in_folder=PRE_KRAKEN2 / "samples",
        out_folder=lambda w: PRE_KRAKEN2 / w.kraken2_db,
        kraken_db_name=lambda w: w.kraken2_db,
        samples=" ".join(SAMPLES),
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2> {log} 1>&2

            mkdir \
                --parents \
                --verbose \
                /dev/shm/{params.kraken_db_name} \
            2>> {log} 1>&2

            mkdir \
                --parents \
                --verbose \
                {params.out_folder} \
            2>> {log} 1>&2

            rsync \
                --archive \
                --progress \
                --recursive \
                --times \
                --verbose \
                --chown $(whoami):$(whoami) \
                --chmod u+rw \
                {input.database}/*.k2d \
                /dev/shm/{params.kraken_db_name} \
            2>> {log} 1>&2

            ( parallel \
                --jobs {threads} \
                --retries 50 \
                kraken2 \
                    --db /dev/shm/{params.kraken_db_name} \
                    --threads 1 \
                    --gzip-compressed \
                    --paired \
                    --output ">(gzip > {params.out_folder}/{{}}.out.gz)" \
                    --report {params.out_folder}/{{}}.k2report \
                    --memory-mapping \
                    {params.in_folder}/{{}}_1.fq.gz \
                    {params.in_folder}/{{}}_2.fq.gz \
                "2>" {params.out_folder}/{{}}.log \
            ::: {params.samples} \
            )

        }} || {{
            echo "Failed job" 2>> {log} 1>&2
            echo "Hostname was $(hostname)" 2>> {log} 1>&2
        }}

        rm \
            --force \
            --recursive \
            --verbose \
            /dev/shm/{params.kraken_db_name} \
        2>>{log} 1>&2

        echo "Finished kraken2 in $(hostname)" 2>> {log} 1>&2
        """


rule preprocess__kraken2__all:
    input:
        out_gzs=[
            PRE_KRAKEN2 / kraken2_db / f"{sample_id}.out.gz"
            for sample_id in SAMPLES
            for kraken2_db in KRAKEN2_DBS
        ],
        reports=[
            PRE_KRAKEN2 / kraken2_db / f"{sample_id}.k2report"
            for sample_id in SAMPLES
            for kraken2_db in KRAKEN2_DBS
        ],
