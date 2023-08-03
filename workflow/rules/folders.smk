READS = Path("results/reads/")

# preprocessing
PRE = Path("results/preprocessing/")
FASTP = PRE / "fastp/"
PRE_INDEX = PRE / "index/"
PRE_BOWTIE2 = PRE / "bowtie2"
NONHOST = PRE / "nonhost/"

# preprocessing evaluation
NONPAREIL = PRE / "nonpareil/"
SINGLEM = PRE / "singlem/"
PRE_COVERM = PRE / "coverm/"


# assembly
ASSEMBLY = Path("results/assembly/")
MEGAHIT = ASSEMBLY / "megahit/"
ASSEMBLY_RENAME = ASSEMBLY / "renaming/"
ASSEMBLY_INDEX = ASSEMBLY / "index/"
ASSEMBLY_BOWTIE2 = ASSEMBLY / "bowtie2/"

# Assembly evaluation
ASSEMBLY_QUAST = ASSEMBLY / "quast/"
ASSEMBLY_COVERM = ASSEMBLY / "coverm/"


# binners
BIN = Path("results/binning/")
VAMB = BIN / "vamb/"  # This could be a metabinner
CONCOCT = BIN / "concoct/"
METABAT2 = BIN / "metabat2/"
MAXBIN2 = BIN / "maxbin2/"
BIN_METAWRAP = BIN / "metawrap/"

# bin evaluation?


# metabinners
METABIN = Path("results/metabinning/")
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
