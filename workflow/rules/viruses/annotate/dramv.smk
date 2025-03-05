rule viruses__annotate__dramv__setup:
    input:
        dram_db=features["databases"]["dram"],
    output:
        setup=touch(VIR_DRAMV / "setup.done"),
    log:
        VIR_DRAMV / "setup.log",
    conda:
        "../../../environments/dram.yml"
    shell:
        """
        DRAM-setup.py set_database_locations \
            --amg_database_loc          {input.dram_db}/amg_database.*.tsv \
            --dbcan_fam_activities_loc  {input.dram_db}/CAZyDB.*.fam-activities.txt \
            --dbcan_loc                 {input.dram_db}/dbCAN-HMMdb-V*.txt \
            --dbcan_subfam_ec_loc       {input.dram_db}/CAZyDB.*.fam.subfam.ec.txt \
            --description_db_loc        {input.dram_db}/description_db.sqlite \
            --etc_module_database_loc   {input.dram_db}/etc_mdoule_database.*.tsv \
            --function_heatmap_form_loc {input.dram_db}/function_heatmap_form.*.tsv \
            --genome_summary_form_loc   {input.dram_db}/genome_summary_form.*.tsv \
            --kofam_hmm_loc             {input.dram_db}/kofam_profiles.hmm \
            --kofam_ko_list_loc         {input.dram_db}/kofam_ko_list.tsv \
            --module_step_form_loc      {input.dram_db}/module_step_form.*.tsv \
            --peptidase_loc             {input.dram_db}/peptidases.*.mmsdb \
            --pfam_hmm_loc              {input.dram_db}/Pfam-A.hmm.dat.gz \
            --pfam_loc                  {input.dram_db}/pfam.mmspro \
            --viral_loc                 {input.dram_db}/refseq_viral.*.mmsdb \
            --vog_annotations_loc       {input.dram_db}/vog_annotations_latest.tsv.gz \
            --vogdb_loc                 {input.dram_db}/vog_latest_hmms.txt \
        2> {log} 1>&2
        """

rule viruses__annotate__dramv__annotate:
    input:
        fa=VIR_VIRSORTER2 / "final-viral-combined-for-dramv.fa.gz",
        tsv=VIR_VIRSORTER2 / "viral-affi-contigs-for-dramv.tab.gz",
        dram_db=features["databases"]["dram"],
        setup=VIR_DRAMV / "setup.done",
    output:
        annotations=VIR_DRAMV / "annotations.tsv.gz",
    log:
        VIR_DRAMV / "annotate.log",
    conda:
        "../../../environments/dram.yml"
    params:
        workdir=VIR_DRAMV,
    threads: 24
    resources:
        mem_mb=32 * 1024,
        time_min=60,
    shell:
        """
        seqtk split \
            -n {threads} \
            {params.workdir}/splits \
            {input.fa} \
        2>> {log} 1>&2

        parallel \
            --jobs {threads} \
            DRAM-v.py annotate \
                --input_fasta {{.}}.fa \
                --output_dir {params.workdir}/{{/.}} \
                --skip_trnascan \
                --virsorter_affi_contigs <(gzip -dc {input.tsv}) \
        ::: {params.workdir}/splits.*.fa \
        2>> {log} 1>&2

        Rscript --vanilla workflow/scripts/stack_dram_annotations.R \
            --output-file {output.annotations} \
            {params.workdir}/splits.*/annotations.tsv \
        2>> {log} 1>&2

        rm \
            --recursive \
            --force \
            --verbose \
            {params.workdir}/splits* \
        2>> {log} 1>&2
        """


rule viruses__annotate__dramv__distill:
    input:
        annotations=VIR_DRAMV / "annotations.tsv.gz",
    output:
        amg_summary=VIR_DRAMV / "amg_summary.tsv.gz",
        vmag_stats=VIR_DRAMV / "vMAG_stats.tsv.gz",
        product=VIR_DRAMV / "product.html",
    log:
        VIR_DRAMV / "distill.log",
    conda:
        "../../../environments/dram.yml"
    params:
        outdir=VIR_DRAMV,
        workdir=VIR_DRAMV / "tmp",
    shadow:
        "minimal"
    shell:
        """
        DRAM-v.py distill \
            --input_file {input.annotations} \
            --output_dir {params.workdir} \
        2> {log} 1>&2

        mv \
            {params.workdir}/* \
            {params.outdir}/ \
        2>> {log} 1>&2

        bgzip \
            --threads {threads} \
            {params.outdir}/amg_summary.tsv \
            {params.outdir}/vMAG_stats.tsv \
        2>> {log} 1>&2
        """


rule viruses__annotate__dramv__all:
    input:
        rules.viruses__annotate__dramv__distill.output,
