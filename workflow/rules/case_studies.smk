# Case study: rank genes inside published linkage intervals (config/loci.tsv).
rule locus_case_study:
    input:
        "config/loci.tsv",
        "results/backtest/interval_simulation.tsv",
    output:
        "results/case_studies/locus_ranking.tsv",
        "results/case_studies/locus_top10.tsv",
    shell:
        "python workflow/scripts/locus_case_study.py"
