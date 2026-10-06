# Stage 1: truth set validation.
rule validate_truth:
    input:
        config["truth"]["path"],
    output:
        "results/tables/truth_validated.txt",
    shell:
        "pnprank validate-truth {input} && echo ok > {output}"
