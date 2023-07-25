READS = Path("results/reads/")

PRE = Path("results/preprocessing/")
FASTP = PRE / "fastp/"
BOWTIE2_PRE = PRE / "bowtie2/"
NONHOST = PRE / "nonhost/"
NONPAREIL = PRE / "nonpareil/"
SINGLEM = PRE / "singlem/"
COVERM = PRE / "coverm/"


ASSEMBLY = Path("results/assembly/")
MEGAHIT = ASSEMBLY / "megahit/"
QUAST = ASSEMBLY / "quast/"
BOWTIE2_INDEXES_ASSEMBLY = ASSEMBLY / "bowtie2_indexes/"
BOWTIE2_ASSEMBLY = ASSEMBLY / "bowtie2/"
