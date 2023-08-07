READS = Path("results/reads/")

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


# binners
BIN = Path("results/bin/")
VAMB = BIN / "vamb/"  # This could be a metabinner
CONCOCT = BIN / "concoct/"
METABAT2 = BIN / "metabat2/"
MAXBIN2 = BIN / "maxbin2/"
BIN_METAWRAP = BIN / "metawrap/"

# bin evaluation?


# metabinners
METABIN = Path("results/metabin/")
MAGSCOT = METABIN / "magscot/"
PRODIGAL = MAGSCOT / "prodigal/"
METABIN_METAWRAP = METABIN / "metawrap/"
METABIN_RENAME = METABIN / "renaming/"
METABIN_INDEX = BIN / "index/"
METABIN_BOWTIE2 = METABIN / "bowtie2/"

# metabinning evaluation
METABIN_QUAST = BIN / "quast/"
METABIN_COVERM = METABIN / "coverm/"
METABIN_GTDBTK = METABIN / "gtdbtk/"
METABIN_DRAM = METABIN / "dram/"


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
