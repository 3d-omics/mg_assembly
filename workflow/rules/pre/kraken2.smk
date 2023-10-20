rule pre_eval_kraken2_assign_all:
    input:
        files=[
            FASTP / f"{sample}.{library}_{ending}.fq.gz"
            for sample, library in SAMPLE_LIBRARY
            for ending in ["1", "2"]
        ],
        database=get_kraken2_database,
    output:
        out_gzs=[
            KRAKEN2 / "{kraken_db}" / f"{sample}.{library}.out.gz"
            for sample, library in SAMPLE_LIBRARY
        ],
        reports=[
            KRAKEN2 / "{kraken_db}" / f"{sample}.{library}.report"
            for sample, library in SAMPLE_LIBRARY
        ],
    log:
        KRAKEN2 / "{kraken_db}.log",
    threads: 24
    resources:
        mem_mb=params["pre"]["kraken2"]["memory_gb"] * 1024,
        runtime=60,
    params:
        in_folder=FASTP,
        out_folder=compose_out_folder_for_eval_kraken2_assign_all,
        kraken_db_shm="/dev/shm/{kraken_db}",
    conda:
        "pre.yml"
    shell:
        """
        mapfile -t sample_ids < <(find "{params.in_folder}" -name "*_1.fq.gz" -exec basename {{}} _1.fq.gz \;)

        {{
            mkdir --parents {params.kraken_db_shm}
            mkdir --parents {params.out_folder}

            rsync \
                -Pravt \
                {input.database}/*.k2d \
                {params.kraken_db_shm} \
            2> {log} 1>&2

            for sample_id in ${{sample_ids[@]}} ; do \

                kraken2 \
                    --db {params.kraken_db_shm} \
                    --threads {threads} \
                    --gzip-compressed \
                    --output >(pigz > {params.out_folder}/${{sample_id}}.out.gz) \
                    --report {params.out_folder}/${{sample_id}}.report \
                    --memory-mapping \
                    {params.in_folder}/${{sample_id}}_1.fq.gz \
                    {params.in_folder}/${{sample_id}}_2.fq.gz \
                2> {params.out_folder}/${{sample_id}}.log  1>&2

            done
        }} || {{
            echo "Failed job" 2>> {log} 1>&2
        }}

        rm -rfv {params.kraken_db_shm} 2>>{log} 1>&2
        """


rule pre_eval_kraken2:
    input:
        [
            KRAKEN2 / f"{kraken_db}/{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken_db in features["kraken2_dbs"]
        ],