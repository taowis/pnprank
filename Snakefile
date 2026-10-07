# pnprank workflow. One rules file per stage under workflow/rules/.
configfile: "config/config.yaml"
configfile: "config/paths.yaml"

include: "workflow/rules/truth.smk"
include: "workflow/rules/features.smk"


rule all:
    input:
        "results/tables/truth_validated.txt",
        "results/tables/stage4/cv_performance.tsv",
        "results/backtest/interval_summary.tsv",
        "results/backtest/repeated_cv_tier1.tsv",
