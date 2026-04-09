rule csvtk__concat:
    input:
        ["filein1.tsv.gz", "filein2.tsv.gz"]
    output:
        "fileout.tsv.gz"
    log:
        "out.log"
    params:
        subcommand="concat",
        extra="--tabs --out-tabs"
    wrapper:
        "v5.2.1/utils/csvtk"


rule csvtk__join__left:
    input:
        ["filein1.tsv.gz", "filein2.tsv.gz"]
    output:
        "fileout.tsv.gz"
    log:
        "log",
    params:
        subcommand="join",
        extra="--left-join --tabs --out-tabs",
    wrapper:
        "v5.2.1/utils/csvtk"
