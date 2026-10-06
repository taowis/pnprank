"""`pnprank` command-line entry point. Subcommands are added as stages complete."""

from __future__ import annotations

import argparse
import sys

from pnprank import __version__
from pnprank.truth import build_candidates, fetch_sources, validate_truth


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pnprank", description=__doc__)
    parser.add_argument("--version", action="version", version=f"pnprank {__version__}")
    sub = parser.add_subparsers(dest="command")

    p_truth = sub.add_parser("validate-truth", help="check the truth TSV schema")
    p_truth.add_argument("path", nargs="?", default="data/truth/neuropathy_genes.tsv")

    p_cand = sub.add_parser("build-candidates", help="snapshot sources and list candidate genes")
    p_cand.add_argument("--config", default="config/config.yaml")
    p_cand.add_argument("--no-download", action="store_true", help="reuse existing snapshots")

    args = parser.parse_args(argv)
    if args.command == "validate-truth":
        problems = validate_truth(args.path)
        for p in problems:
            print(p, file=sys.stderr)
        return 1 if problems else 0
    if args.command == "build-candidates":
        return _build_candidates(args)
    parser.print_help()
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
