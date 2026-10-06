# pnprank workflow. One rules file per stage under workflow/rules/.
configfile: "config/config.yaml"
configfile: "config/paths.yaml"

include: "workflow/rules/truth.smk"


rule all:
    input:
        "results/tables/truth_validated.txt",
