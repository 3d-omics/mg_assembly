rule virify_download_db:
    output:
        db_folder=features["virify_database"],
    log:
        VIRIFY / "databases.log",
    conda:
        "virify.yml"
    threads: 1
    shell:
        """
        nextflow run EBI-Metagenomics/emg-viral-pipeline \
            -revision v2.0.0 \
            --databases {output.db_folder} \
            --fasta /dev/null \
            --cores {threads} \
            --max_cores {threads} \
            --workdir {output.db_folder}/work \
            -profile local,singularity \
        2> {log} 1>&2

        rm -rf results/null
        """


rule virify_analyze_one:
    input:
        fasta=ASSEMBLE_RENAME / "{assembly_id}.fa",
        db_folder=features["virify_database"],
    output:
        out_folder=directory(VIRIFY / "{assembly_id}"),
    log:
        VIRIFY / "{assembly_id}.log",
    conda:
        "virify.yml"
    threads: 4
    params:
        mem_gb=12,
        output=VIRIFY,
    resources:
        runtime=6 * 60,
        mem_mb=12 * 1024,
    shell:
        """
        nextflow run EBI-Metagenomics/emg-viral-pipeline \
            -revision v2.0.0 \
            --fasta {input.fasta} \
            --output {params.output} \
            --cores {threads} \
            --max_cores {threads} \
            --memory {params.mem_gb} \
            --workdir {output.out_folder}/work \
            --databases {input.db_folder} \
            --singularity_cachedir {output.out_folder}/singularity \
            -profile local,docker \
        2> {log} 1>&2
        """


rule virify_all:
    input:
        [VIRIFY / f"{assembly_id}" for assembly_id in ASSEMBLIES],
    conda:
        "virify.yml"
    log:
        VIRIFY / "all.log",
    shell:
        """
        nextflow clean -f 2> {log} 1>&2
        """


rule virify:
    input:
        rules.virify_all.input,


localrules:
    virify_download_db,
    virify_one,
