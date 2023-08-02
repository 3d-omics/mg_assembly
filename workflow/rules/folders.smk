READS = Path("results/reads/")

PRE = Path("results/preprocessing/")
FASTP = PRE / "fastp/"
BOWTIE2_PRE = PRE / "bowtie2/"
NONHOST = PRE / "nonhost/"
NONPAREIL = PRE / "nonpareil/"
SINGLEM = PRE / "singlem/"
COVERM_PRE = PRE / "coverm/"


ASSEMBLY = Path("results/assembly/")
MEGAHIT = ASSEMBLY / "megahit/"
MEGAHIT_RENAMING = ASSEMBLY / "renaming/"
QUAST = ASSEMBLY / "quast/"
BOWTIE2_INDEXES_ASSEMBLY = ASSEMBLY / "indexes/"
BOWTIE2_ASSEMBLY = ASSEMBLY / "bowtie2/"
COVERM_ASSEMBLY = ASSEMBLY / "coverm/"

BINNING = Path("results/binning/")
# binners
VAMB = BINNING / "vamb/"
CONCOCT = BINNING / "concoct/"
METABAT2 = BINNING / "metabat2/"
MAXBIN2 = BINNING / "maxbin2/"
METAWRAP_BINNING = BINNING / "binning/"


# metabinners
MAGSCOT = BINNING / "magscot/"
MAGSCOT_PRODIGAL = MAGSCOT / "prodigal/"
METAWRAP_REFINEMENT = BINNING / "refinement/"
METAWRAP_RENAMING = BINNING / "renaming/"

# binning evaluation
BINNING_QUAST = BINNING / "quast/"
BOWTIE2_INDEXES_BINNING = BINNING / "indexes/"
BOWTIE2_BINNING = BINNING / "bowtie2/"
COVERM_BINNING = BINNING / "coverm/"
GTDBTK = BINNING / "gtdbtk/"
DRAM = BINNING / "dram/"
DRAM_ANNOTATE = DRAM / "annotate/"
DRAM_DISTILL = DRAM / "distill/"

# dereplicate
DEREPLICATE = Path("results/dereplicate/")
DREP = DEREPLICATE / "drep/"
