# rule _assemble__concoct__cut_up_fasta:
#     """Run concoct cut_up_fasta.py on one assembly."""
#     input:
#         assembly=MEGAHIT / "{assembly_id}.fa.gz",
#     output:
#         assembly_10k=CONCOCT / "{assembly_id}_data" / "cut.fa.gz",
#         bed_10k=CONCOCT / "{assembly_id}_data" / "cut.bed",
#     log:
#         CONCOCT / "{assembly_id}_data" / "cut.log",
#     conda:
#         "concoct.yml"
#     shell:
#         """
#         ( cut_up_fasta.py \
#             <(gzip --decompress --stdout {input.assembly}) \
#             --chunk_size 10000 \
#             --overlap_size 0 \
#             --merge_last \
#             --bedfile {output.bed_10k} \
#         | gzip \
#             --best \
#         > {output.assembly_10k} \
#         ) 2> {log}
#         """


# rule _assemble__concoct__coverage_table:
#     """Run concoct concoct_coverage_table.py on one assembly."""
#     input:
#         crams=get_crams_from_assembly_id,
#         crais=get_crais_from_assembly_id,
#         assembly=MEGAHIT / "{assembly_id}.fa.gz",
#         bed_10k=CONCOCT / "{assembly_id}_data" / "cut.bed",
#     output:
#         coverage=CONCOCT / "{assembly_id}_data" / "coverage.tsv.gz",
#     log:
#         CONCOCT / "{assembly_id}_data" / "coverage.log",
#     conda:
#         "concoct.yml"
#     resources:
#         mem_mb=8 * 1024,
#     params:
#         bams=compose_bams_for_concoct_coverage_table,
#         bais=compose_bais_for_concoct_coverage_table,
#     shell:
#         """
#         for cram in {input.crams} ; do

#             bam={CONCOCT}/$(basename $cram .cram).bam

#             samtools view \
#                 --exclude-flags 4 \
#                 --fast \
#                 --output $bam \
#                 --output-fmt BAM \
#                 --reference {input.assembly} \
#                 --threads {threads} \
#                 $cram

#             samtools index $bam

#         done 2> {log} 1>&2

#         concoct_coverage_table.py \
#             {input.bed_10k} \
#             {params.bams} \
#         | gzip -9 \
#         > {output.coverage} \
#         2>> {log}

#         rm \
#             --verbose \
#             --force \
#             {params.bams} \
#             {params.bais} \
#         2>> {log} 1>&2
#         """


# rule _assemble__concoct__run:
#     """Run concoct on one assembly."""
#     input:
#         assembly_10k=CONCOCT / "{assembly_id}_data" / "cut.fa.gz",
#         coverage=CONCOCT / "{assembly_id}_data" / "coverage.tsv.gz",
#     output:
#         args=CONCOCT / "{assembly_id}_data" / "run_args.txt",
#         clustering=CONCOCT / "{assembly_id}_data" / "run_clustering_gt1000.csv",
#         log=CONCOCT / "{assembly_id}_data" / "run_log.txt",
#         original_data=CONCOCT / "{assembly_id}_data" / "run_original_data_gt1000.csv",
#         components=CONCOCT / "{assembly_id}_data" / "run_PCA_components_data_gt1000.csv",
#         transformed_data=CONCOCT
#         / "{assembly_id}_data"
#         / "run_PCA_transformed_data_gt1000.csv",
#     log:
#         CONCOCT / "{assembly_id}_data" / "run.log",
#     conda:
#         "concoct.yml"
#     params:
#         basename=lambda w: CONCOCT / f"{w.assembly_id}_data/run",
#     threads: 24
#     resources:
#         runtime=24 * 60,
#         mem_mb=double_ram(4),
#     retries: 5
#     shell:
#         """
#         concoct \
#             --threads {threads} \
#             --composition_file <(gzip -dc {input.assembly_10k}) \
#             --coverage_file <(gzip -dc {input.coverage}) \
#             --basename {params.basename} \
#         2>> {log} 1>&2
#         """


# rule _assemble__concoct__merge_cutup_clustering:
#     """Run concoct merge_cutup_clustering.py on one assembly."""
#     input:
#         clustering=CONCOCT / "{assembly_id}_data" / "run_clustering_gt1000.csv",
#     output:
#         clustering_merged=CONCOCT / "{assembly_id}_data" / "merge.csv",
#     log:
#         CONCOCT / "{assembly_id}_data" / "merge.log",
#     conda:
#         "concoct.yml"
#     shell:
#         """
#         merge_cutup_clustering.py \
#             {input.clustering} \
#         > {output.clustering_merged} \
#         2>> {log}
#         """


# rule _assemble__concoct__extract_fasta_bins:
#     """Run concoct extract_fasta_bins.py on one assembly."""
#     input:
#         assembly=MEGAHIT / "{assembly_id}.fa.gz",
#         clustering_merged=CONCOCT / "{assembly_id}_data" / "merge.csv",
#     output:
#         bins=directory(CONCOCT / "{assembly_id}_bins"),
#     log:
#         CONCOCT / "{assembly_id}_bins.log",
#     conda:
#         "concoct.yml"
#     resources:
#         mem_mb=4 * 1024,
#     shell:
#         """
#         mkdir --parents {output.bins} --verbose 2> {log} 1>&2

#         extract_fasta_bins.py \
#             <(gzip --decompress --stdout {input.assembly}) \
#             {input.clustering_merged} \
#             --output_path {output.bins} \
#         2>> {log} 1>&2

#         find \
#             {output.bins} \
#             -name "*.fa" \
#             -exec pigz --best --verbose {{}} \; \
#         2>> {log} 1>&2
#         """


rule _assemble__concoct:
    input:
        assembly=MEGAHIT / "{assembly_id}.fa.gz",
        crams=get_crams_from_assembly_id,
    output:
        directory(CONCOCT / "{assembly_id}")
    log:
        CONCOCT / "{assembly_id}.log"
    conda:
        "concoct.yml"
    resources:
        mem_mb=double_ram(8)
    # retries:
    #     5
    params:
        workdir=lambda w: CONCOCT / w.assembly_id
    shell:
        """
        mkdir --parents --verbose {params.workdir} 2> {log} 1>&2

        cut_up_fasta.py \
            <(gzip --decompress --stdout {input.assembly}) \
            --chunk_size 10000 \
            --overlap_size 0 \
            --merge_last \
            --bedfile {params.workdir}/cut.bed \
        > {params.workdir}/cut.fa \
        2>> {log}

        for cram in {input.crams} ; do

            bam={params.workdir}/$(basename $cram .cram).bam

            samtools view \
                --exclude-flags 4 \
                --fast \
                --output $bam \
                --output-fmt BAM \
                --reference {input.assembly} \
                --threads {threads} \
                $cram

            samtools index $bam

        done 2>> {log} 1>&2

        concoct_coverage_table.py \
            {params.workdir}/cut.bed \
            {params.workdir}/*.bam \
        > {params.workdir}/coverage.tsv \
        2>> {log}

        concoct \
            --threads {threads} \
            --composition_file {params.workdir}/cut.fa \
            --coverage_file {params.workdir}/coverage.tsv \
            --basename {params.workdir}/run \
        2>> {log} 1>&2

        merge_cutup_clustering.py \
            {params.workdir}/run_clustering_gt1000.csv \
        > {params.workdir}/merge.csv \
        2>> {log}

        extract_fasta_bins.py \
            <(gzip --decompress --stdout {input.assembly}) \
            {params.workdir}/merge.csv \
            --output_path {params.workdir} \
        2>> {log} 1>&2

        rm \
            --force \
            --verbose \
            {params.workdir}/cut.fa \
            {params.workdir}/cut.bed \
            {params.workdir}/coverage.tsv \
            {params.workdir}/*.csv \
            {params.workdir}/*.bam \
            {params.workdir}/*.bai \
            {params.workdir}/*.txt \
        2>> {log} 1>&2

        pigz \
            --best \
            --verbose \
            {params.workdir}/*.fa \
        2>> {log} 1>&2
        """


rule assemble__concoct:
    """Run concoct on all assemblies"""
    input:
        [CONCOCT / f"{assembly_id}" for assembly_id in ASSEMBLIES],
