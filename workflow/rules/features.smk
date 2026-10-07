# Stage 4: standardise atlases (raw downloads must already be in place via
# `pnprank download-atlas`; RDS objects are converted by the R scripts first), then
# build features and evaluate models by nested CV.
INTERIM = config["data_root"] + "/interim"


rule standardize_atlas:
    output:
        INTERIM + "/{atlas}/pseudobulk_sums.tsv.gz",
        INTERIM + "/{atlas}/pseudobulk_groups.tsv",
    shell:
        "python workflow/scripts/standardize.py {wildcards.atlas}"


rule stage4_model:
    input:
        expand(INTERIM + "/{atlas}/pseudobulk_groups.tsv", atlas=config["features"]["atlases"]),
        "config/celltype_map.tsv",
    output:
        "results/tables/stage4/cv_performance.tsv",
        "results/tables/stage4/paired_auroc_differences.tsv",
        "results/tables/stage4/ranking_current_logistic_full.tsv.gz",
    shell:
        "python workflow/scripts/stage4_model.py"


# Stage 5: temporal back-test, interval simulation, repeated CV and summaries.
rule stage5_backtest:
    input:
        "results/tables/stage4/cv_performance.tsv",
    output:
        "results/backtest/backtest_performance.tsv",
        "results/backtest/interval_simulation.tsv",
    shell:
        "python workflow/scripts/stage5_backtest.py 5"


rule stage5_repeated_cv:
    input:
        "results/tables/stage4/cv_performance.tsv",
    output:
        "results/backtest/repeated_cv_tier1.tsv",
    shell:
        "python workflow/scripts/stage5_repeated_cv.py 4"


rule stage5_summarise:
    input:
        "results/backtest/backtest_performance.tsv",
        "results/backtest/interval_simulation.tsv",
    output:
        "results/backtest/interval_summary.tsv",
        "results/backtest/failure_analysis.tsv",
    shell:
        "python workflow/scripts/stage5_summarise.py"
