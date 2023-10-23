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

# Assembly evaluation
ASSEMBLE_QUAST = ASSEMBLE / "quast/"
ASSEMBLE_COVERM = ASSEMBLE / "coverm/"


# Viral prediction
VIRIFY = Path("results/virify/")

# binners
BIN = Path("results/bin/")
VAMB = BIN / "vamb/"  # This could be a metabinner
CONCOCT = BIN / "concoct/"
METABAT2 = BIN / "metabat2/"
MAXBIN2 = BIN / "maxbin2/"
MAGSCOT = BIN / "magscot/"
PRODIGAL = MAGSCOT / "prodigal/"
BIN_QUAST = BIN / "quast/"


# dereplicate
DEREPLICATE = Path("results/dereplicate/")
DREP = DEREPLICATE / "drep/"
DREP_INDEX = DEREPLICATE / "index/"
DREP_BOWTIE2 = DEREPLICATE / "bowtie2/"


# dereplicate evaluation
DREP_GTDBTK = DEREPLICATE / "gtdbtk/"
DREP_QUAST = DEREPLICATE / "quast/"
DREP_COVERM = DEREPLICATE / "coverm/"
DREP_DRAM = DEREPLICATE / "dram/"


# reports
REPORT = Path("report/")
REPORT_STEP = REPORT / "by_step/"
REPORT_ASSEMBLY = REPORT / "by_assembly/"
