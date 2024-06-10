READS = Path("results/reads/")

# reference
REFERENCE = Path("results/reference/")
HOSTS = REFERENCE / "hosts"

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

# Prokaryotes
PROK = Path("results/prokaryotes")

PROK_CLUSTER = PROK / "cluster"
ASSEMBLE_INDEX = PROK_CLUSTER / "index/"
ASSEMBLE_BOWTIE2 = PROK_CLUSTER / "bowtie2/"
CONCOCT = PROK_CLUSTER / "concoct/"
METABAT2 = PROK_CLUSTER / "metabat2/"
MAXBIN2 = PROK_CLUSTER / "maxbin2/"
MAGSCOT = PROK_CLUSTER / "magscot/"
PRODIGAL = PROK_CLUSTER / "prodigal/"
DREP = PROK_CLUSTER / "drep/"


# quantify
PROK_QUANT = PROK / "quantify"
QUANT_INDEX = PROK_QUANT / "index/"
QUANT_BOWTIE2 = PROK_QUANT / "bowtie2/"
COVERM = PROK_QUANT / "coverm/"


# dereplicate evaluation
PROK_ANN = PROK / "annotate"
GTDBTK = PROK_ANN / "gtdbtk/"
QUAST = PROK_ANN / "quast/"
DRAM = PROK_ANN / "dram/"
CHECKM = PROK_ANN / "checkm2"


# Viral prediction
VIR = Path("results/viruses/")

VIR_CLUSTER = VIR / "cluster"
GENOMADC = VIR_CLUSTER / "genomad"
CHECKVC = VIR_CLUSTER / "checkv"
DEDUPE = VIR_CLUSTER / "dedupe"
MMSEQS = VIR_CLUSTER / "mmseqs"

VIR_ANN = VIR / "annotate"
VIRSORTER2 = VIR_ANN / "virsorter2"
GENOMADA = VIR_ANN / "genomad"
DRAMV = VIR_ANN / "dramv"
QUASTV = VIR_ANN / "quast"
CHECKV = VIR_ANN / "checkv"

VIR_QUANT = VIR / "quantify"
VINDEX = VIR_QUANT / "index"
VBOWTIE2 = VIR_QUANT / "bowtie2"
VCOVERM = VIR_QUANT / "coverm"


# reports
REPORT = Path("reports/")
REPORT_STEP = REPORT / "step/"
REPORT_SAMPLE = REPORT / "library/"
