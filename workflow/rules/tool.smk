# Stage 6: the ranking shipped with the `pnprank` CLI (same out-of-fold procedure as the
# Stage 5 back-test, applied to today's tier-1 labels).
rule stage6_ranking:
    input:
        "results/tables/stage4/cv_performance.tsv",
        config["truth"]["path"],
    output:
        "results/tables/pnprank_ranking.tsv.gz",
        "results/tables/pnprank_ranking_meta.json",
    shell:
        "pnprank build-ranking --n-jobs 5 --out {output[0]}"


rule stage6_checks:
    input:
        "results/tables/pnprank_ranking.tsv.gz",
    output:
        "results/tables/stage6/release_checks.tsv",
    shell:
        "python workflow/scripts/stage6_checks.py"


rule gene_coordinates:
    input:
        "results/tables/pnprank_ranking.tsv.gz",
        ancient(config["external"]["human_gtf"]["dest"]),
    output:
        config["external"]["human_gtf"]["coordinates_table"],
    shell:
        "python workflow/scripts/gene_coordinates.py"
