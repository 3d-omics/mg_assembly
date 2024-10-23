include: "dram_functions.smk"

rule prokaryotes__annotate__dram__setup:
    """
    Set up the databases from DRAM, no matter what the config file says.
    """
    input:
        dram_db=features["databases"]["dram"],
    output:
        touch(PROK_ANN / "dram.setup.txt"),
    log:
        PROK_ANN / "dram.setup.log",
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
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__dram__annotate:
    """Annotate dereplicate genomes with DRAM"""
    input:
        fasta=MAGS / "{mag_id}.fa",
        dram_db=features["databases"]["dram"],
        setup=PROK_ANN / "dram.setup.txt",
    output:
        work_dir=temp(directory(PROK_ANN / "dram.annotate" / "{mag_id}")),
    log:
        PROK_ANN / "dram.annotate" / "{mag_id}.log",
    conda:
        "../../../environments/dram.yml"
    params:
        min_contig_size=params["prokaryotes"]["annotate"]["dram"]["annotate"][
            "min_contig_size"
        ],
    shell:
        """
        rm \
            --recursive \
            --force \
            --verbose \
            {output.work_dir} \
        2> {log} 1>&2

        DRAM.py annotate \
            --input_fasta {input.fasta} \
            --output_dir {output.work_dir} \
            --threads 1 \
        2>> {log} 1>&2
        """




rule prokaryotes__annotate__dram__annotate__aggregate_annotations:
    """Aggregate DRAM annotations"""
    input:
        collect_dram_annotate,
    output:
        PROK_ANN / "dram.annotations.tsv.gz",
    log:
        PROK_ANN / "dram.annotate.aggregate.log",
    conda:
        "../../../environments/dram.yml"
    params:
        work_dir=PROK_ANN / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/annotations.tsv \
        | csvformat \
            --out-tabs \
        | bgzip \
            --compress-level 9 \
        > {output} ) \
        2> {log}
        """


rule prokaryotes__annotate__dram__annotate__aggregate_trnas:
    """Aggregate DRAM tRNAs"""
    input:
        collect_dram_annotate,
    output:
        PROK_ANN / "dram.trnas.tsv",
    log:
        PROK_ANN / "dram.trnas.log",
    conda:
        "../../../environments/dram.yml"
    params:
        work_dir=PROK_ANN / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/trnas.tsv \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule prokaryotes__annotate__dram__annotate_aggregate_rrnas:
    """Aggregate DRAM rRNAs"""
    input:
        collect_dram_annotate,
    output:
        PROK_ANN / "dram.rrnas.tsv",
    log:
        PROK_ANN / "dram.rrnas.log",
    conda:
        "../../../environments/dram.yml"
    params:
        work_dir=PROK_ANN / "dram.annotate",
    shell:
        """
        ( csvstack \
            --tabs \
            {params.work_dir}/*/rrnas.tsv \
        | csvformat \
            --out-tabs \
        > {output} ) \
        2> {log}
        """


rule prokaryotes__annotate__dram__annotate_archive:
    """
    Create tarball once annotations are merged done
    """
    input:
        work_dirs=collect_dram_annotate,
        annotations=PROK_ANN / "dram.annotations.tsv.gz",
        trnas=PROK_ANN / "dram.trnas.tsv",
        rrnas=PROK_ANN / "dram.rrnas.tsv",
    output:
        tarball=PROK_ANN / "dram.annotate.tar.gz",
    log:
        PROK_ANN / "dram.archive.log",
    conda:
        "../../../environments/dram.yml"
    params:
        out_dir=PROK_ANN,
        work_dir=PROK_ANN / "dram.annotate",
    shell:
        """
        tar \
            --create \
            --file {output.tarball} \
            --use-compress-program="pigz --processes {threads}" \
            --verbose \
            {params.work_dir} \
        2>> {log} 1>&2

        rm -rfv {params.work_dir}
        """


rule prokaryotes__annotate__dram__distill:
    """Distill DRAM annotations."""
    input:
        annotations=PROK_ANN / "dram.annotations.tsv.gz",
        trnas=PROK_ANN / "dram.trnas.tsv",
        rrnas=PROK_ANN / "dram.rrnas.tsv",
        dram_db=features["databases"]["dram"],
        setup=PROK_ANN / "dram.setup.txt",
    output:
        work_dir=temp(directory(PROK_ANN / "dram.distill")),
    log:
        PROK_ANN / "dram.distill.log",
    conda:
        "../../../environments/dram.yml"
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --output_dir {output.work_dir} \
            --rrna_path  {input.rrnas} \
            --trna_path  {input.trnas} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__dram__distill_archive:
    input:
        work_dir=PROK_ANN / "dram.distill",
    output:
        genome=PROK_ANN / "dram.genome_stats.tsv",
        metabolism=PROK_ANN / "dram.metabolism_summary.xlsx",
        product_tsv=PROK_ANN / "dram.product.tsv",
    log:
        PROK_ANN / "dram.distill_archive.log",
    conda:
        "../../../environments/dram.yml"
    params:
        out_dir=PROK_ANN,
    shell:
        """
        for file in genome_stats.tsv metabolism_summary.xlsx product.tsv ; do

            cp \
                --verbose \
                {input.work_dir}/$file \
                {params.out_dir}/dram.$file \

        done 2> {log} 1>&2
        """


rule prokaryotes__annotate__dram__all:
    """Run DRAM on dereplicated genomes."""
    input:
        rules.prokaryotes__annotate__dram__annotate_archive.output,
        rules.prokaryotes__annotate__dram__distill_archive.output,


localrules:
    prokaryotes__annotate__dram__distill_archive,
