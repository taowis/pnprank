"""`pnprank` command-line entry point. Subcommands are added as stages complete."""

from __future__ import annotations

import argparse
import sys

from pnprank import __version__
from pnprank.truth import validate_truth


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pnprank", description=__doc__)
    parser.add_argument("--version", action="version", version=f"pnprank {__version__}")
    sub = parser.add_subparsers(dest="command")

    p_truth = sub.add_parser("validate-truth", help="check the truth TSV schema")
    p_truth.add_argument("path", nargs="?", default="data/truth/neuropathy_genes.tsv")

    args = parser.parse_args(argv)
    if args.command == "validate-truth":
        problems = validate_truth(args.path)
        for p in problems:
            print(p, file=sys.stderr)
        return 1 if problems else 0
    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
