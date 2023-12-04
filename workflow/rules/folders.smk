READS = Path("results/reads/")
REFERENCE = Path("results/reference/")

# preprocess
PRE = Path("results/preprocess/")
FASTP = PRE / "fastp/"
PRE_INDEX = PRE / "index/"
PRE_BOWTIE2 = PRE / "bowtie2"
NONHOST = PRE / "nonhost/"

# preprocess evaluation
NONPAREIL = PRE / "nonpareil/"
SINGLEM = PRE / "singlem/"
PRE_COVERM = PRE / "coverm/"
KRAKEN2 = PRE / "kraken2/"


# assemble
ASSEMBLE = Path("results/assemble/")
MEGAHIT = ASSEMBLE / "megahit/"
ASSEMBLE_RENAME = ASSEMBLE / "renaming/"
ASSEMBLE_INDEX = ASSEMBLE / "index/"
ASSEMBLE_BOWTIE2 = ASSEMBLE / "bowtie2/"
CONCOCT = ASSEMBLE / "concoct/"
METABAT2 = ASSEMBLE / "metabat2/"
MAXBIN2 = ASSEMBLE / "maxbin2/"
MAGSCOT = ASSEMBLE / "magscot/"
PRODIGAL = MAGSCOT / "prodigal/"
DREP = ASSEMBLE / "drep/"

# VAMB = ASSEMBLE / "vamb/"  # This could be a metabinner


# Viral prediction
VIRIFY = Path("results/virify/")


# dereplicate - this is confusing
DEREPLICATE = Path("results/dereplicate/")
DREP_INDEX = DEREPLICATE / "index/"
DREP_BOWTIE2 = DEREPLICATE / "bowtie2/"


# dereplicate evaluation
DREP_GTDBTK = DEREPLICATE / "gtdbtk/"
DREP_QUAST = DEREPLICATE / "quast/"
DREP_COVERM = DEREPLICATE / "coverm/"
DREP_DRAM = DEREPLICATE / "dram/"
DREP_CHECKM = DEREPLICATE / "checkm2"

# reports
REPORT = Path("report/")
REPORT_STEP = REPORT / "by_step/"
REPORT_ASSEMBLY = REPORT / "by_assembly/"
