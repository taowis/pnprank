# Comparators inside the temporal back-test (docs/stage_reports/comparators.md):
# GO guilt-by-association and BioGRID random walk, each on releases dated before the
# cutoff (config/comparators.yaml), plus gnomAD constraint scores.
CMP_EXT = config["data_root"] + "/external/comparators"
CMP_OUT = "results/backtest/comparators"


# ancient(): editing comments in the config should not trigger a re-download check.
rule comparators_fetch:
    input:
        ancient("config/comparators.yaml"),
    output:
        CMP_EXT + "/manifest.json",
    shell:
        "python workflow/scripts/comparators_fetch.py"


# Reads the per-cutoff PNPrank scores cached by stage5_backtest.
rule comparators_backtest:
    input:
        CMP_EXT + "/manifest.json",
        "results/backtest/backtest_performance.tsv",
        ancient("logs/leakage_tests_passed.txt"),
    output:
        CMP_OUT + "/comparator_performance.tsv",
        CMP_OUT + "/comparator_paired_differences.tsv",
        CMP_OUT + "/comparator_coverage.tsv",
        CMP_OUT + "/comparator_data_coverage.tsv",
        CMP_OUT + "/comparator_random_seed_null.tsv",
    shell:
        "python workflow/scripts/comparators_backtest.py results/backtest/scores 5"


rule comparators_summarise:
    input:
        CMP_OUT + "/comparator_performance.tsv",
        CMP_OUT + "/comparator_paired_differences.tsv",
    output:
        CMP_OUT + "/comparators_auroc.png",
        CMP_OUT + "/comparator_headline_differences.tsv",
    shell:
        "python workflow/scripts/comparators_summarise.py"
