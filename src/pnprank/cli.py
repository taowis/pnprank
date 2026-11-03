"""`pnprank` command-line entry point. Subcommands are added as stages complete."""

from __future__ import annotations

import argparse
import sys

from pathlib import Path

from pnprank import __version__
from pnprank.truth import assemble_truth, build_candidates, fetch_sources, validate_truth

DEFAULT_RANKING = Path(__file__).resolve().parents[2] / "results" / "tables" / "pnprank_ranking.tsv.gz"
DEFAULT_COORDS = DEFAULT_RANKING.with_name("gene_coordinates_e116.tsv.gz")
DISPLAY = ["rank", "hgnc_symbol", "status", "percentile", "percentile_seed_min",
           "percentile_seed_max", "percentile_expr_matched", "percentile_non_atlas",
           "gtex_mean_decile", "top_cell_type"]
NOTE = ("note: percentiles rank genes for follow-up; they are not evidence that a gene causes "
        "neuropathy. Known tier-1 genes were training positives and have no percentile.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pnprank", description=__doc__)
    parser.add_argument("--version", action="version", version=f"pnprank {__version__}")
    sub = parser.add_subparsers(dest="command")

    p_truth = sub.add_parser("validate-truth", help="check the truth TSV schema")
    p_truth.add_argument("path", nargs="?", default="data/truth/neuropathy_genes.tsv")

    p_cand = sub.add_parser("build-candidates", help="snapshot sources and list candidate genes")
    p_cand.add_argument("--config", default="config/config.yaml")
    p_cand.add_argument("--no-download", action="store_true", help="reuse existing snapshots")

    p_asm = sub.add_parser("assemble-truth", help="check curated papers on PubMed, write truth TSV")
    p_asm.add_argument("--config", default="config/config.yaml")
    p_asm.add_argument("--curation", default="data/truth/curation.json")
    p_asm.add_argument("--review", default="data/truth/review.tsv")

    p_dl = sub.add_parser("download-atlas", help="download listed files of one atlas")
    p_dl.add_argument("atlas_id")
    p_dl.add_argument("files", nargs="+")
    p_dl.add_argument("--base-url", required=True)
    p_dl.add_argument("--registry", default="config/atlases.yaml")
    p_dl.add_argument("--paths", default="config/paths.yaml")

    p_build = sub.add_parser("build-ranking", help="fit the release model and write the ranking")
    p_build.add_argument("--config", default="config/config.yaml")
    p_build.add_argument("--paths", default="config/paths.yaml")
    p_build.add_argument("--out", default=str(DEFAULT_RANKING))
    p_build.add_argument("--n-jobs", type=int, default=1)

    p_q = sub.add_parser("query", help="look up genes in the ranking")
    p_q.add_argument("genes", nargs="*", help="HGNC symbols (aliases accepted) or Ensembl IDs")
    p_q.add_argument("--file", help="file with one gene per line ('-' for stdin)")
    p_q.add_argument("--ranking", default=str(DEFAULT_RANKING))
    p_q.add_argument("--tsv", action="store_true", help="tab-separated output, all columns")

    p_top = sub.add_parser("top", help="list the highest-ranked genes")
    p_top.add_argument("-n", type=int, default=50)
    p_top.add_argument("--include-known", action="store_true",
                       help="include known tier-1/2 neuropathy genes")
    p_top.add_argument("--expr-matched", action="store_true",
                       help="order by the expression-matched percentile")
    p_top.add_argument("--ranking", default=str(DEFAULT_RANKING))
    p_top.add_argument("--tsv", action="store_true", help="tab-separated output, all columns")

    p_int = sub.add_parser(
        "interval", help="rank the genes inside one or more GRCh38 intervals",
        description="Rank the protein-coding genes overlapping each interval. Coordinates are "
                    "GRCh38, 1-based and inclusive (e.g. chr7:150,000,000-160,000,000, "
                    "optionally named: DHMN1=chr7:...). Lift intervals on other assemblies "
                    "over to GRCh38 first.")
    p_int.add_argument("regions", nargs="*", help="[name=]chrom:start-end")
    p_int.add_argument("--bed", help="BED file of intervals (0-based starts, as BED defines)")
    p_int.add_argument("--build", default="GRCh38", help="assembly of the coordinates "
                       "(only GRCh38/hg38 is accepted)")
    p_int.add_argument("--pad-kb", type=float, default=0,
                       help="extend each interval by this many kb on both sides")
    p_int.add_argument("--expr-matched", action="store_true",
                       help="order by the expression-matched percentile")
    p_int.add_argument("--ranking", default=str(DEFAULT_RANKING))
    p_int.add_argument("--coordinates", default=str(DEFAULT_COORDS))
    p_int.add_argument("--tsv", action="store_true", help="tab-separated output, all columns")

    args = parser.parse_args(argv)
    if args.command == "validate-truth":
        problems = validate_truth(args.path)
        for p in problems:
            print(p, file=sys.stderr)
        return 1 if problems else 0
    if args.command == "build-candidates":
        return _build_candidates(args)
    if args.command == "assemble-truth":
        return _assemble_truth(args)
    if args.command == "download-atlas":
        return _download_atlas(args)
    if args.command == "build-ranking":
        return _build_ranking(args)
    if args.command in ("query", "top"):
        return _lookup(args)
    if args.command == "interval":
        return _interval(args)
    parser.print_help()
    return 0


def _download_atlas(args: argparse.Namespace) -> int:
    from pathlib import Path

    import yaml

    from pnprank.atlases import download_atlas_files

    paths = yaml.safe_load(Path(args.paths).read_text())
    for rec in download_atlas_files(
        args.registry, args.atlas_id, args.files, args.base_url,
        paths["data_root"], paths["disk_limit_gb"],
    ):
        print(f"{rec['name']}\t{rec['size_bytes']}\t{rec['sha256']}")
    return 0


def _write_tsv(path, rows, fieldnames) -> None:
    import csv

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _assemble_truth(args: argparse.Namespace) -> int:
    import csv
    import json
    from pathlib import Path

    import yaml

    from pnprank.truth import REQUIRED_COLUMNS

    cfg = yaml.safe_load(Path(args.config).read_text())["truth"]
    with open(cfg["candidates"], newline="") as fh:
        candidates = list(csv.DictReader(fh, delimiter="\t"))
    curation = json.loads(Path(args.curation).read_text())
    truth, review = assemble_truth(curation, candidates, frozenset(cfg.get("force_tier2", [])))
    _write_tsv(cfg["path"], truth, REQUIRED_COLUMNS)
    _write_tsv(args.review, review, ["hgnc_symbol", "flag", "detail"])
    print(f"{len(truth)} truth rows -> {cfg['path']}; {len(review)} review flags -> {args.review}")
    return 0


def _build_candidates(args: argparse.Namespace) -> int:
    import csv
    from pathlib import Path

    import yaml

    cfg = yaml.safe_load(Path(args.config).read_text())["truth"]
    panels = cfg["panelapp_panels"]
    panel_ids = {"hsp" if k == "hsp" else k: v for k, v in panels.items()}
    srcdir = Path(cfg["sources_dir"])
    if not args.no_download:
        fetch_sources(srcdir, panel_ids)
    rows = build_candidates(srcdir, panel_ids)
    out = Path(cfg["candidates"])
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} candidates -> {out}")
    return 0


def _build_ranking(args: argparse.Namespace) -> int:
    import json

    import pandas as pd
    import yaml

    from pnprank.features import label_features, load_celltype_map, load_profiles
    from pnprank.ranking import build_ranking

    cfg = yaml.safe_load(Path(args.config).read_text())
    root = Path(yaml.safe_load(Path(args.paths).read_text())["data_root"])
    fc, rc = cfg["features"], cfg["release"]
    X = pd.read_csv(root / "processed" / "features_static.tsv.gz", sep="\t", index_col=0)
    universe = list(X.index)
    profiles = load_profiles(fc["atlases"], root / "interim", load_celltype_map(fc["celltype_map"]),
                              cfg["pseudobulk"]["min_cells_per_group"],
                              cfg["pseudobulk"]["min_donors"],
                         root / "processed" / "features_static.tsv.gz")
    truth = pd.read_csv(cfg["truth"]["path"], sep="\t", dtype=str)
    cands = pd.read_csv(cfg["truth"]["candidates"], sep="\t", dtype=str)
    hgnc = pd.read_csv(Path(cfg["truth"]["sources_dir"]) / "hgnc_complete_set.txt", sep="\t",
                       dtype=str)
    hgnc = hgnc[hgnc.ensembl_gene_id.isin(universe) & (hgnc.status == "Approved")]
    symbols = dict(zip(hgnc.ensembl_gene_id, hgnc.symbol))
    aliases = {g: "|".join(filter(None, [a or "", p or ""])) for g, a, p in
               zip(hgnc.ensembl_gene_id, hgnc.alias_symbol.fillna(""), hgnc.prev_symbol.fillna(""))}
    seeds = [cfg["seeds"]["global"] + i for i in range(rc["stability_seeds"])]
    tab, meta = build_ranking(
        X, truth, set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id), symbols, aliases,
        lambda known: label_features(profiles, universe, known),
        [c for c in X.columns if c.startswith(("gnomad|", "gtex|"))],
        X[rc["expression_covariate"]], seeds, rc["scoring_folds"], cfg["model"]["inner_folds"],
        rc["expression_bins"], args.n_jobs)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out, sep="\t", index=False, compression="gzip")
    meta = {"pnprank_version": __version__, "n_genes": len(tab), **meta}
    meta_path = out.with_name(out.name.removesuffix(".tsv.gz") + "_meta.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"{len(tab)} genes -> {out}; metadata -> {meta_path}")
    return 0


def _lookup(args: argparse.Namespace) -> int:
    from pnprank.ranking import load_ranking, resolve, top_genes

    path = Path(args.ranking)
    if not path.exists():
        print(f"ranking not found: {path} (run `pnprank build-ranking` first)", file=sys.stderr)
        return 2
    tab = load_ranking(path)
    if args.command == "top":
        out = top_genes(tab, args.n, args.include_known,
                        "percentile_expr_matched" if args.expr_matched else "percentile")
        problems = {}
    else:
        queries = list(args.genes)
        if args.file:
            fh = sys.stdin if args.file == "-" else open(args.file)
            queries += [ln.split("\t")[0].strip() for ln in fh if not ln.startswith("#")]
        if not queries:
            print("no genes given", file=sys.stderr)
            return 2
        out, problems = resolve(tab, queries)
    if args.tsv:
        out.to_csv(sys.stdout, sep="\t", index=False)
    else:
        cols = (["query"] if "query" in out.columns else []) + DISPLAY
        print(out[cols].to_string(index=False, na_rep="-") if len(out) else "(no matches)")
    for q, why in problems.items():
        print(f"{q}: {why}", file=sys.stderr)
    print(NOTE, file=sys.stderr)
    return 1 if problems and out.empty else 0


def _interval(args: argparse.Namespace) -> int:
    import pandas as pd

    from pnprank.ranking import (GRCH38_NAMES, genes_in_interval, load_coordinates,
                                 load_ranking, parse_region, read_bed)

    if args.build.lower() not in GRCH38_NAMES:
        print(f"build {args.build!r} not supported: the ranking uses GRCh38 (Ensembl 116) gene "
              "coordinates. Lift the intervals over to GRCh38 (e.g. UCSC liftOver) and check "
              "them against the flanking markers.", file=sys.stderr)
        return 2
    try:
        regions = [parse_region(r) for r in args.regions]
        if args.bed:
            with open(args.bed) as fh:
                regions += read_bed(fh)
    except (ValueError, IndexError) as err:
        print(err, file=sys.stderr)
        return 2
    if not regions:
        print("no intervals given", file=sys.stderr)
        return 2
    for p in (args.ranking, args.coordinates):
        if not Path(p).exists():
            print(f"not found: {p}", file=sys.stderr)
            return 2
    tab, coords = load_ranking(args.ranking), load_coordinates(args.coordinates)
    by = "percentile_expr_matched" if args.expr_matched else "percentile"
    pad = int(args.pad_kb * 1000)
    blocks = []
    for name, chrom, start, end in regions:
        out = genes_in_interval(tab, coords, chrom, start, end, pad, by)
        out.insert(0, "interval", name)
        blocks.append(out)
        if args.tsv:
            continue
        n_known = int(out.interval_rank.isna().sum())
        span = f"chr{chrom}:{max(1, start - pad):,}-{end + pad:,}"
        print(f"== {name}  {span} (GRCh38)  {len(out)} genes, {len(out) - n_known} scored, "
              f"{n_known} known tier-1")
        cols = ["interval_rank"] + DISPLAY + ["start"]
        print(out[cols].to_string(index=False, na_rep="-") if len(out) else "(no genes)")
        print()
    if args.tsv:
        pd.concat(blocks, ignore_index=True).to_csv(sys.stdout, sep="\t", index=False)
    print(NOTE + " interval_rank orders scored genes within the interval (1 = highest); "
          "genes partly overlapping an edge are included.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
