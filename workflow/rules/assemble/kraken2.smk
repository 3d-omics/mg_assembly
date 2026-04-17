rule assemble__kraken2__assign_contigs:
    input:
        fastas=[ASMB_MEGAHIT / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
    output:
        out_gzs=[
            ASMB_KRAKEN2 / "{kraken2_db}" / f"{assembly_id}.out.gz"
            for assembly_id in ASSEMBLIES
        ],
        reports=[
            ASMB_KRAKEN2 / "{kraken2_db}" / f"{assembly_id}.report"
            for assembly_id in ASSEMBLIES
        ],
    log:
        ASMB_KRAKEN2 / "{kraken2_db}.log",
    conda:
        ENVS / "kraken2.yml"
    threads: 8
    resources:
        mem_mb=800 * 1024,
        runtime=6 * 60,
    params:
        in_folder=ASMB_MEGAHIT,
        out_folder=lambda w: ASMB_KRAKEN2 / w.kraken2_db,
        kraken_db_name=lambda w: w.kraken2_db,
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2>> {log} 1>&2

            mkdir --parents /dev/shm/{params.kraken_db_name}
            mkdir --parents {params.out_folder}

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

            for file in {input.fastas} ; do \

                sample_id=$(basename $file .fa.gz)
                fasta={params.in_folder}/${{sample_id}}.fa.gz
                output={params.out_folder}/${{sample_id}}.out.gz
                report={params.out_folder}/${{sample_id}}.report
                log={params.out_folder}/${{sample_id}}.log

                echo $(date) Processing $sample_id 2>> {log} 1>&2

                kraken2 \
                    --db /dev/shm/{params.kraken_db_name} \
                    --threads {threads} \
                    --gzip-compressed \
                    --output >(pigz --processes {threads} --best > $output) \
                    --report $report \
                    --memory-mapping \
                    $fasta \
                2> $log 1>&2

            done
        }} || {{
            echo "Failed job" 2>> {log} 1>&2
        }}

        rm \
            --force \
            --recursive \
            --verbose \
            /dev/shm/{params.kraken_db_name} \
        2>>{log} 1>&2
        """


rule assemble__kraken2__assign_contigs__all:
    input:
        [
            ASMB_KRAKEN2 / f"{kraken2_db}" / f"{assembly_id}.out.gz"
            for assembly_id in ASSEMBLIES
            for kraken2_db in features["databases"]["kraken2"]
        ],
        [
            ASMB_KRAKEN2 / f"{kraken2_db}" / f"{assembly_id}.report"
            for assembly_id in ASSEMBLIES
            for kraken2_db in features["databases"]["kraken2"]
        ],


rule assemble__kraken2__all:
    input:
        rules.assemble__kraken2__assign_contigs__all.input,
