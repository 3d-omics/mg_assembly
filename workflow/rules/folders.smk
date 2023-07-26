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
QUAST = ASSEMBLY / "quast/"
BOWTIE2_INDEXES_ASSEMBLY = ASSEMBLY / "bowtie2_assembly_indexes/"
BOWTIE2_ASSEMBLY = ASSEMBLY / "bowtie2_assembly/"
COVERM_ASSEMBLY = ASSEMBLY / "coverm_assembly/"

BINNING = Path("results/binning/")
METAWRAP_BINNING = BINNING / "metawrap_binning/"
METAWRAP_REFINEMENT = BINNING / "metawrap_refinement/"
METAWRAP_RENAMING = BINNING / "metawrap_renaming/"
BOWTIE2_INDEXES_BINNING = BINNING / "bowtie2_bins_indexes/"
BOWTIE2_BINNING = BINNING / "bowtie2_bins/"
COVERM_BINNING = BINNING / "coverm_bins/"
