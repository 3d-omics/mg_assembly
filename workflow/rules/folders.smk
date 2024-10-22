# reference
# REFERENCE = Path("results/reference/")
# HOSTS = REFERENCE / "hosts"

RESULTS = Path("results/")

# Preprocess
PRE = RESULTS / "preprocess"
PRE_FASTP = PRE / "fastp/"
PRE_INDEX = PRE / "index/"
PRE_BOWTIE2 = PRE / "bowtie2"

# Assemble
ASSEMBLE = RESULTS / "assemble"
ASSEMBLE_MEGAHIT = ASSEMBLE / "megahit/"
ASSEMBLE_INDEX = ASSEMBLE / "index/"
ASSEMBLE_BOWTIE2 = ASSEMBLE / "bowtie2/"
ASSEMBLE_QUAST = ASSEMBLE / "quast/"

# Prokaryotes
PROK = RESULTS / "prokaryotes"

## Prokaryotes - Cluster
PROK_CLUSTER = PROK / "cluster"
CONCOCT = PROK_CLUSTER / "concoct/"
METABAT2 = PROK_CLUSTER / "metabat2/"
MAXBIN2 = PROK_CLUSTER / "maxbin2/"
MAGSCOT = PROK_CLUSTER / "magscot/"
PRODIGAL = PROK_CLUSTER / "prodigal/"
DREP = PROK_CLUSTER / "drep/"


## Prokaryotes - Quantify
PROK_QUANT = PROK / "quantify"
QUANT_INDEX = PROK_QUANT / "index/"
QUANT_BOWTIE2 = PROK_QUANT / "bowtie2/"
COVERM = PROK_QUANT / "coverm/"

## Prokaryotes - Annotate
PROK_ANN = PROK / "annotate"
MAGS = PROK_ANN / "mags/"
GTDBTK = PROK_ANN / "gtdbtk/"
QUAST = PROK_ANN / "quast/"
DRAM = PROK_ANN / "dram/"


# Viruses
VIR = RESULTS / "viruses"

## Viruses - Cluster
VIR_CLUSTER = VIR / "cluster"
GENOMADC = VIR_CLUSTER / "genomad"
CHECKVC = VIR_CLUSTER / "checkv"
DEDUPE = VIR_CLUSTER / "dedupe"
MMSEQS = VIR_CLUSTER / "mmseqs"

## Viruses - Annotation
VIR_ANN = VIR / "annotate"
VIRSORTER2 = VIR_ANN / "virsorter2"
GENOMADA = VIR_ANN / "genomad"
DRAMV = VIR_ANN / "dramv"
QUASTV = VIR_ANN / "quast"
CHECKV = VIR_ANN / "checkv"

## Viruses - Quantify
VIR_QUANT = VIR / "quantify"
VINDEX = VIR_QUANT / "index"
VBOWTIE2 = VIR_QUANT / "bowtie2"
VCOVERM = VIR_QUANT / "coverm"


# reports
REPORT = Path("reports/")
REPORT_STEP = REPORT / "step/"
REPORT_SAMPLE = REPORT / "library/"
