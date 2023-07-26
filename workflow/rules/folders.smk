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
METAWRAP_BINNING = ASSEMBLY / "metawrap_binning/"
METAWRAP_REFINEMENT = ASSEMBLY / "metawrap_refinement/"
METAWRAP_RENAMING = ASSEMBLY / "metawrap_renaming/"
BOWTIE2_INDEXES_BINS = ASSEMBLY / "bowtie2_bins_indexes/"
BOWTIE2_BINS = ASSEMBLY / "bowtie2_bins/"
COVERM_ASSEMBLY = ASSEMBLY / "coverm_assembly/"
COVERM_BINS = ASSEMBLY / "coverm_bins/"
