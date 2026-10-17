# pnprank workflow. Rules are grouped by stage under workflow/rules/; docs/workflow.md
# describes the order, inputs and outputs of every step.
configfile: "config/config.yaml"
configfile: "config/paths.yaml"

include: "workflow/rules/external.smk"
include: "workflow/rules/truth.smk"
include: "workflow/rules/features.smk"
include: "workflow/rules/comparators.smk"
include: "workflow/rules/tool.smk"
include: "workflow/rules/case_studies.smk"


# Everything that does not depend on live web queries (see case_studies.smk).
rule all:
    default_target: True
    input:
        "results/tables/truth_validated.txt",
        "results/pilot/summary.json",
        "results/tables/stage4/cv_performance.tsv",
        "results/backtest/interval_summary.tsv",
        "results/backtest/repeated_cv_tier1.tsv",
        "results/backtest/comparators/comparator_headline_differences.tsv",
        "results/tables/stage6/release_checks.tsv",
        "results/tables/gene_coordinates_e116.tsv.gz",
        "results/case_studies/locus_ranking.tsv",
        "results/case_studies/patient_tissue_expression.tsv",
