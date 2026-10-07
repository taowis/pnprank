# pnprank workflow. One rules file per stage under workflow/rules/.
configfile: "config/config.yaml"
configfile: "config/paths.yaml"

include: "workflow/rules/truth.smk"
include: "workflow/rules/features.smk"


rule all:
    input:
        "results/tables/truth_validated.txt",
        "results/tables/stage4/cv_performance.tsv",
