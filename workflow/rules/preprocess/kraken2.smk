rule preprocess__kraken2__assign:
    """
    Run kraken2 over all samples at once using the /dev/shm/ trick.

    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=[
            PRE_FASTP / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        rerverses=[
            PRE_FASTP / f"{sample_id}.{library_id}_2.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
    output:
        out_gzs=[
            PRE_KRAKEN2 / "{kraken2_db}" / f"{sample_id}.{library_id}.out.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        reports=[
            PRE_KRAKEN2 / "{kraken2_db}" / f"{sample_id}.{library_id}.k2report"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    log:
        PRE_KRAKEN2 / "{kraken2_db}.log",
    params:
        in_folder=PRE_FASTP,
        out_folder=lambda w: PRE_KRAKEN2 / w.kraken2_db,
        kraken_db_name=lambda w: w.kraken2_db,
        sample_libs=[
            f"{sample_id}.{library_id}" for sample_id, library_id in SAMPLE_LIBRARY
        ],
    threads: 24
    resources:
        mem_mb=2 * 800 * 1024,  # Use twice the size of the database, we use /dev/shm
        runtime=24 * 60,
    conda:
        "../../environments/kraken2.yml"
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2> {log} 1>&2

            mkdir \
                --parents \
                --verbose \
                /dev/shm/{params.kraken_db_name} \
            2>> {log} 1>&2

            mkdir --parents --verbose {params.out_folder} 2>> {log} 1>&2

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
            ::: {params.sample_libs} \
            )

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


rule preprocess__kraken2__bracken:
    input:
        database=lambda w: features["databases"]["kraken2"][w.kraken2_db],
        report=PRE_KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.k2report",
    output:
        bracken=touch(
            PRE_KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.{level}.bracken"
        ),
    log:
        PRE_KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.{level}.bracken.log",
    conda:
        "../../environments/kraken2.yml"
    params:
        extra=params["preprocess"]["kraken2"]["bracken"]["extra"],
        level=lambda w: w.level,
    shell:
        """
        if [ ! -s {input.report} ] ; then
            echo "Empty report. Skipping" 2> {log} 1>&2
            exit 0
        fi

        bracken \
            -d {input.database} \
            -i {input.report} \
            -o {output.bracken} \
            -l {params.level} \
            {params.extra} \
        2> {log} 1>&2
        """


rule preprocess__kraken2__bracken__combine:
    """Combine all the bracken outputs for a single database"""
    input:
        lambda w: [
            PRE_KRAKEN2 / w.kraken2_db / f"{sample_id}.{library_id}.{w.level}.bracken"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        PRE_KRAKEN2 / "{kraken2_db}.{level}.tsv",
    log:
        PRE_KRAKEN2 / "{kraken2_db}.{level}tsv.log",
    conda:
        "../../environments/kraken2.yml"
    shell:
        """
        combine_bracken_outputs.py \
            --files {input} \
            --output {output} \
        2> {log} 1>&2
        """


rule preprocess__kraken2__all:
    """Get the combined bracken results for all databases"""
    input:
        [
            PRE_KRAKEN2 / f"{kraken2_db}.{level}.tsv"
            for kraken2_db in KRAKEN2_DBS
            for level in "DPCOFGS".split("")
        ],
