# Stages 3-5. Atlases are standardised (RDS objects are converted to Matrix Market by
# the R scripts, called from standardize.py) and pseudobulked, then features are built
# and models evaluated by nested CV, then back-tested.
INTERIM = config["data_root"] + "/interim"
HGNC = config["truth"]["sources_dir"] + "/hgnc_complete_set.txt"


# Reference tables are fixed by release and checksummed (fetch_external), so they are
# marked ancient(): rebuilding them identically does not force every atlas to rerun.
rule standardize_atlas:
    input:
        lambda w: atlas_raw_files(w.atlas),
        ancient(HGNC),
        ancient(config["external"]["mouse_gtf"]["genes_table"]),
        ancient(config["reference"]["orthology"]["table"]),
    output:
        INTERIM + "/{atlas}/pseudobulk_sums.tsv.gz",
        INTERIM + "/{atlas}/pseudobulk_groups.tsv",
    shell:
        "python workflow/scripts/standardize.py {wildcards.atlas}"


# Stage 3: feasibility pilot on one atlas (all discovery years; not a back-test).
# Reads the count matrix that standardize_atlas writes for this atlas (counts/).
rule stage3_pilot:
    input:
        INTERIM + "/" + config["pilot"]["atlas"] + "/counts/matrix.mtx.gz",
        ancient(INTERIM + "/" + config["pilot"]["atlas"] + "/pseudobulk_groups.tsv"),
    output:
        "results/pilot/summary.json",
        "results/pilot/auroc_by_celltype.tsv",
    shell:
        "python workflow/scripts/pilot.py"


# Stage 4: full feature matrix and nested-CV model comparison on current labels.
rule stage4_model:
    input:
        expand(INTERIM + "/{atlas}/pseudobulk_groups.tsv", atlas=config["features"]["atlases"]),
        "config/celltype_map.tsv",
        config["external"]["gnomad_constraint"]["dest"],
        config["external"]["gtex_median_tpm"]["dest"],
        config["truth"]["path"],
    output:
        "results/tables/stage4/cv_performance.tsv",
        "results/tables/stage4/paired_auroc_differences.tsv",
        "results/tables/stage4/ranking_current_logistic_full.tsv.gz",
    shell:
        "python workflow/scripts/stage4_model.py"


# Stage 5: temporal back-test, interval simulation, repeated CV and summaries.
# stage5_backtest also caches per-cutoff out-of-fold scores in results/backtest/scores/
# (not committed), which the comparator and locus rules read.
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


# Leakage control: tests/test_leakage.py must pass before back-test results are
# summarised. ancient() so that re-running the tests does not force a re-summary.
rule leakage_tests:
    input:
        "tests/test_leakage.py",
    output:
        "logs/leakage_tests_passed.txt",
    shell:
        "python -m pytest -q tests/test_leakage.py > {output}"


rule stage5_summarise:
    input:
        "results/backtest/backtest_performance.tsv",
        "results/backtest/interval_simulation.tsv",
        ancient("logs/leakage_tests_passed.txt"),
    output:
        "results/backtest/interval_summary.tsv",
        "results/backtest/failure_analysis.tsv",
    shell:
        "python workflow/scripts/stage5_summarise.py"
