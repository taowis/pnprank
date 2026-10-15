# Case study: rank genes inside published linkage intervals (config/loci.tsv).
# Interval gene lists are fetched from the Ensembl REST API on first run and cached
# under <data_root>/external/loci/.
rule locus_case_study:
    input:
        "config/loci.tsv",
        "results/backtest/interval_simulation.tsv",
        config["external"]["ucsc_cytoband"]["dest"],
    output:
        "results/case_studies/locus_ranking.tsv",
        "results/case_studies/locus_top10.tsv",
    shell:
        "python workflow/scripts/locus_case_study.py"


rule patient_tissue_expression:
    input:
        "data/truth/neuropathy_genes.tsv",
        ancient(config["external"]["gtex_median_tpm"]["dest"]),
    output:
        "results/case_studies/patient_tissue_expression.tsv",
        "results/case_studies/patient_tissue_expression_genes.tsv",
    shell:
        "python workflow/scripts/patient_tissue_expression.py"


# The two rules below query live web services (PubMed, HGNC, PanelApp, ClinGen,
# Ensembl, Alliance). Results change as those services are updated, so they are not
# part of `all`; run them by name (`snakemake -c1 live_queries`). Responses for the
# C. elegans table are cached under <data_root>/external/celegans/ and reused.
rule top50_literature_queries:
    input:
        "results/tables/pnprank_ranking.tsv.gz",
    output:
        "results/tables/top50_literature_queries.tsv",
        "results/tables/top50_literature_curated.tsv",
    shell:
        "python workflow/scripts/top50_literature_queries.py -n 50"


rule celegans_orthologues:
    input:
        "results/tables/pnprank_ranking.tsv.gz",
    output:
        "results/tables/celegans_shortlist.tsv",
        "results/tables/tier1_celegans_orthologues.tsv",
    shell:
        "python workflow/scripts/celegans_orthologues.py"


rule live_queries:
    input:
        rules.top50_literature_queries.output,
        rules.celegans_orthologues.output,
