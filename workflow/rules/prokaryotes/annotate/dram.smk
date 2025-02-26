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
        fasta=PROK_MAGS / "{mag_id}.fa",
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
    resources:
        mem_mb=8 * 1024,
        runtime=24 * 60,
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


for file in ["annotations", "trnas", "rrnas"]:

    rule:
        name:
            f"prokaryotes__annotate__dram__annotate__aggregate_{file}"
        input:
            collect_dram_annotate,
        output:
            PROK_ANN / f"dram.{file}.tsv.gz",
        log:
            PROK_ANN / f"dram.{file}.log",
        conda:
            "../../../environments/dram.yml"
        params:
            work_dir=PROK_ANN / "dram.annotate",
        threads: 24
        shell:
            f"( csvtk concat --tabs {{params.work_dir}}/*/{file} "
            f"| sed -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' "
            f"| bgzip --compress-level 9 --threads {{threads}} "
            f"> {{output}} "
            f") 2> {{log}}"


for file in ["genes.gff", "genes.fna", "genes.faa", "scaffolds.fna"]:

    rule:
        name:
            f"prokaryotes__annotate__dram__annotate__concatenate_{file}"
        input:
            collect_dram_annotate,
        output:
            PROK_ANN / f"dram.{file}.gz",
        log:
            PROK_ANN / f"dram.{file}.log",
        conda:
            "../../../environments/dram.yml"
        params:
            work_dir=PROK_ANN / "dram.annotate",
        threads: 24
        shell:
            f"( cat {{params.work_dir}}/*/{file} "
            f"| sed -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' "
            f"| bgzip --compress-level 9 --threads {{threads}}"
            f"> {{output}}"
            f") 2> {{log}}"


rule prokaryotes__annotate__dram__annotate__aggregate_genbank:
    """Aggregate all DRAM genbank files"""
    input:
        collect_dram_annotate,
    output:
        PROK_ANN / "dram.genbank.gbk.gz",
    log:
        PROK_ANN / "dram.genbank.log",
    conda:
        "../../../environments/dram.yml"
    params:
        work_dir=PROK_ANN / "dram.annotate",
    threads: 24
    shell:
        """
        ( cat \
            {params.work_dir}/*/genbank/*.gbk \
        | sed \
            -r 's/[[:alnum:]]+:bin_[0-9]+_([[:alnum:]]+:bin_[0-9]+@contig_[0-9]+)/\1/g' \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
        > {output} \
        ) 2> {log}
        """


rule prokaryotes__annotate__dram__annotate__archive:
    """
    Create tarball once annotations are merged done
    """
    input:
        annotations=PROK_ANN / "dram.annotations.tsv.gz",
        trnas=PROK_ANN / "dram.trnas.tsv.gz",
        rrnas=PROK_ANN / "dram.rrnas.tsv.gz",
        gtf=PROK_ANN / "dram.genes.gff.gz",
        fna=PROK_ANN / "dram.genes.fna.gz",
        faa=PROK_ANN / "dram.genes.faa.gz",
        scaffolds=PROK_ANN / "dram.scaffolds.fna.gz",
        genbank=PROK_ANN / "dram.genbank.gbk.gz",
    output:
        tarball=PROK_ANN / "dram.annotate.tar.gz",
    log:
        PROK_ANN / "dram.archive.log",
    conda:
        "../../../environments/dram.yml"
    threads: 24
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
        trnas=PROK_ANN / "dram.trnas.tsv.gz",
        rrnas=PROK_ANN / "dram.rrnas.tsv.gz",
        dram_db=features["databases"]["dram"],
        setup=PROK_ANN / "dram.setup.txt",
    output:
        work_dir=temp(directory(PROK_ANN / "dram.distill")),
    log:
        PROK_ANN / "dram.distill.log",
    conda:
        "../../../environments/dram.yml"
    resources:
        mem_mb=16 * 1024,
        runtime=24 * 60,
    shell:
        """
        DRAM.py distill \
            --input_file {input.annotations} \
            --output_dir {output.work_dir} \
            --rrna_path  {input.rrnas} \
            --trna_path  {input.trnas} \
        2>> {log} 1>&2
        """


rule prokaryotes__annotate__dram__distill__archive:
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
    threads: 24
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
        rules.prokaryotes__annotate__dram__annotate__archive.output,
        rules.prokaryotes__annotate__dram__distill__archive.output,


localrules:
    prokaryotes__annotate__dram__distill__archive,
