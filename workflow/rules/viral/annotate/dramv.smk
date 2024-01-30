rule _viral__annotate__dramv:
    input:
        fa=VIRSORTER2 / "for-dramv" / "final-viral-combined-for-dramv.fa",
        tsv=VIRSORTER2 / "for-dramv" / "viral-affi-contigs-for-dramv.tab",
        dram_db=features["databases"]["dram"],
    output:
        genome=DRAMV / "genome_stats.tsv",
        metabolism=DRAMV / "metabolism_summary.xlsx",
        product_html=DRAMV / "product.html",
        product_tsv=DRAMV / "product.tsv",
    log:
        DRAMV / "dramv.log",
    conda:
        "__environment__.yml"
    params:
        workdir=DRAMV,
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
        2>> {log} 1>&2

        DRAM-v.py annotate \
            --input_fasta {input.fa} \
            --output_dir {params.workdir} \
            --skip_trnascan \
            --threads {threads} \
            --virsorter_affi_contigs {input.tsv} \
        2>> {log}

        DRAM-v.py distill \
            --input_file {params.workdir}/annotations.tsv \
            --output_dir {params.workdir} \
        2>> {log}
        """


rule viral__annotate__dramv:
    input:
        rules._viral__annotate__dramv.input,
