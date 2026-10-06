"""`pnprank` command-line entry point. Subcommands are added as stages complete."""

from __future__ import annotations

import argparse
import sys

from pnprank import __version__
from pnprank.truth import assemble_truth, build_candidates, fetch_sources, validate_truth


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


if __name__ == "__main__":
    sys.exit(main())
