"""Command-line entry point."""

from __future__ import annotations

import argparse
from pathlib import Path

from tnbc_biomarkers.config import AnalysisConfig
from tnbc_biomarkers.pipeline import run, write_ci_subset


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="TNBC versus non-TNBC transcriptomic analysis of a GEO SuperSeries."
    )
    parser.add_argument("--accession", default="GSE76275", help="GEO series accession to fetch")
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path("data/raw"),
        help="Directory for downloaded GEO files",
    )
    parser.add_argument("--outdir", type=Path, default=Path("results"), help="Output directory")
    parser.add_argument(
        "--write-ci-subset",
        type=Path,
        default=None,
        help="Also write a small cached matrix for tests",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = AnalysisConfig(accession=args.accession)
    summary = run(args.accession, args.cache, args.outdir, config)
    if args.write_ci_subset is not None:
        write_ci_subset(summary, args.write_ci_subset)
        print(f"Wrote {args.write_ci_subset}", flush=True)


if __name__ == "__main__":
    main()
