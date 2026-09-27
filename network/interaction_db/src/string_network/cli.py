"""argparse entry point: ``python -m string_network <command>``."""
from __future__ import annotations

import argparse
import sys

from .build import BuildOptions, build_network, write_output
from .compare import compare_files
from .config import apply_overrides, load_config, resolve_output_path
from .download import download_all


def _add_build_overrides(p: argparse.ArgumentParser) -> None:
    p.add_argument("--threshold", type=int, default=None)
    p.add_argument("--threshold-when", dest="threshold_when",
                   choices=["before", "after"], default=None)
    p.add_argument("--agg", choices=["mean", "max", "min", "median", "first"], default=None)
    p.add_argument("--decimals", type=int, default=None)
    p.add_argument("--mapping-source", dest="mapping_source", default=None)
    p.add_argument("--reviewed-only", dest="reviewed_only", action="store_true", default=None)
    p.add_argument("--keep-list", dest="keep_list", default=None)
    p.add_argument("--links-source", dest="links_source",
                   choices=["full", "physical"], default=None)
    p.add_argument("--score-mode", dest="score_mode",
                   choices=["raw", "physical_else_scaled"], default=None)
    p.add_argument("--scale-factor", dest="scale_factor", type=float, default=None)
    p.add_argument("--mapping-fallback", dest="mapping_fallback",
                   action=argparse.BooleanOptionalAction, default=None,
                   help="Enable/disable the UniProt->preferred_name->raw fallback chain.")
    p.add_argument("--mapping-alias-sources", dest="mapping_alias_sources", default=None,
                   help="Comma-separated STRING alias sources for UniProt coverage, "
                        "e.g. UniProt_AC,Ensembl_UniProt.")
    p.add_argument("--output", default=None, help="Override output parquet path.")


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--config", default="config.yaml")

    ap = argparse.ArgumentParser(prog="string_network", parents=[common])
    sub = ap.add_subparsers(dest="command", required=True)

    d = sub.add_parser("download", parents=[common], help="Download source files.")
    d.add_argument("--force", action="store_true")

    b = sub.add_parser("build", parents=[common], help="Build the output parquet.")
    _add_build_overrides(b)

    c = sub.add_parser("compare", parents=[common], help="Compare build against the reference.")
    c.add_argument("--build", default=None, help="Path to a build parquet.")
    c.add_argument("--threshold", type=int, default=None,
                   help="Build threshold (resolves the build path; guards the report).")

    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_config(args.config)

    if args.command == "download":
        download_all(cfg, force=args.force)
        return 0

    if args.command == "build":
        cfg = apply_overrides(cfg, args)
        opts = BuildOptions.from_config(cfg)
        df, _stats = build_network(cfg, opts)
        out = resolve_output_path(cfg, args.output)
        write_output(df, out)
        return 0

    if args.command == "compare":
        cfg = apply_overrides(cfg, args)
        threshold = int(cfg.get("threshold", 700))
        ref_threshold = 700  # the reference file is a >= 700 network
        build_path = args.build or resolve_output_path(cfg)
        ref_path = cfg["paths"]["reference"]

        if threshold == ref_threshold:
            out_path = f"{cfg['paths']['reports']}/comparison_vs_reference.md"
            title = f"Comparison vs reference (build threshold {threshold})"
            metrics, report = compare_files(build_path, ref_path, out_path, title=title)
            print(report)
            return 0

        # Build threshold != reference threshold: never overwrite the canonical
        # >=700 report. Compare only the build's >=700 subset, clearly labelled.
        import pandas as pd
        print(f"[compare] build threshold is {threshold} but the reference is a "
              f">= {ref_threshold} network; comparing only the build's "
              f">= {ref_threshold} subset (the >=700 report is left untouched).")
        build_df = pd.read_parquet(build_path)
        subset = build_df[build_df["combined_score"] >= ref_threshold].reset_index(drop=True)
        ref_df = pd.read_parquet(ref_path)
        from .compare import compare, format_report
        metrics = compare(subset, ref_df)
        title = (f"Comparison vs reference (build threshold {threshold}; "
                 f">= {ref_threshold} subset only)")
        report = format_report(metrics, title=title)
        out_path = (f"{cfg['paths']['reports']}/"
                    f"comparison_ge{threshold}_subset_vs_reference.md")
        import os
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as fh:
            fh.write(report)
        print(f"[compare] wrote {out_path}")
        print(report)
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
