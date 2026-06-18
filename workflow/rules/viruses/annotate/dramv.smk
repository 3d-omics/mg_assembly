rule viruses__annotate__dramv__setup:
    input:
        dram_db=features["databases"]["dram"],
    output:
        setup=touch(VIR_DRAMV / "setup.done"),
    log:
        VIR_DRAMV / "setup.log",
    conda:
        ENVS / "dram.yml"
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
    """Annotate MVP's representative viral catalog with DRAM-v.

    Fed by MVP_06's --DRAM output (a fasta + an affi-contigs-style tsv mimicking
    VirSorter2's format), which is cohort-wide, so unlike before this is a single
    run rather than one per assembly_id.
    """
    input:
        fasta=get_dram_input_fasta,
        tsv=get_dram_input_tsv,
        dram_db=features["databases"]["dram"],
        setup=VIR_DRAMV / "setup.done",
    output:
        annotations=VIR_DRAMV / "annotate" / "annotations.tsv.gz",
        genes_faa=VIR_DRAMV / "annotate" / "genes.faa.gz",
        genes_fna=VIR_DRAMV / "annotate" / "genes.fna.gz",
        genes_gff=VIR_DRAMV / "annotate" / "genes.gff.gz",
        scaffolds_fna=VIR_DRAMV / "annotate" / "scaffolds.fna.gz",
    log:
        VIR_DRAMV / "annotate.log",
    conda:
        ENVS / "dram.yml"
    threads: 8
    resources:
        mem_mb=double_ram(32 * 1024),
        runtime=24 * 60,
    params:
        workdir=VIR_DRAMV / "annotate",
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {params.workdir} \
        2> {log} 1>&2

        DRAM-v.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {params.workdir} \
            --skip_trnascan \
            --virsorter_affi_contigs {input.tsv} \
            --threads {threads} \
        2>> {log} 1>&2

        bgzip \
            --threads {threads} \
            {params.workdir}/annotations.tsv \
            {params.workdir}/genes.faa \
            {params.workdir}/genes.fna \
            {params.workdir}/genes.gff \
            {params.workdir}/scaffolds.fna \
        2>> {log} 1>&2
        """


rule viruses__annotate__dramv__distill:
    input:
        annotations=VIR_DRAMV / "annotate" / "annotations.tsv.gz",
    output:
        amg_summary=VIR_DRAMV / "amg_summary.tsv.gz",
        vmag_stats=VIR_DRAMV / "vMAG_stats.tsv.gz",
        product=VIR_DRAMV / "product.html",
    log:
        VIR_DRAMV / "distill.log",
    conda:
        ENVS / "dram.yml"
    params:
        outdir=VIR_DRAMV,
        workdir=VIR_DRAMV / "tmp",
    resources:
        mem_mb=8 * 1024,
    shell:
        """
        rm -rfv {params.workdir} 2> {log}

        DRAM-v.py distill \
            --input_file {input.annotations} \
            --output_dir {params.workdir} \
        2>> {log} 1>&2

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
        rules.viruses__annotate__dramv__annotate.output,
        rules.viruses__annotate__dramv__distill.output,
