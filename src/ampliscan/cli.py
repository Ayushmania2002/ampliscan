"""Command-line interface for ampliscan.

Uses only the standard library (argparse) so the CLI works without any extra
dependencies. Subcommands:

    ampliscan gui       -- launch the desktop GUI (no terminal needed)
    ampliscan demux     -- demultiplex FASTQ/FASTA into per-bin files
    ampliscan validate  -- lint a panel YAML and print its summary
    ampliscan sniff     -- preview the bin distribution on the first N reads
    ampliscan panel     -- build a panel YAML from the PI's Excel sheet

Running ``ampliscan`` with no arguments launches the GUI.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from ampliscan import __version__


def _load_panel(config_path: str):
    from ampliscan.config import Panel
    return Panel.from_yaml(config_path)


def _make_policy(args):
    from ampliscan.config import MatchPolicy
    return MatchPolicy(
        anchor_mismatch=args.anchor_mismatch,
        barcode_mismatch=args.barcode_mismatch,
        allow_indels=args.indels,
    )


def cmd_demux(args) -> int:
    from ampliscan import demux_to_bins

    panel = _load_panel(args.config)
    policy = _make_policy(args)
    workers = args.workers  # 0 => auto-detect all cores; default 1

    from ampliscan.io import HAVE_ISAL
    n_cores = (os.cpu_count() or 1) if workers == 0 else workers
    gz = "isal (fast)" if HAVE_ISAL else "stdlib (pip install ampliscan[fast] to speed up)"
    print(f"ampliscan {__version__}  |  panel: {panel.name}  |  cores: {n_cores}  |  gzip: {gz}",
          file=sys.stderr)

    last = [0.0]

    def progress(n):
        now = time.time()
        if now - last[0] > 1.0:
            print(f"\r  processed {n:,} reads...", end="", file=sys.stderr, flush=True)
            last[0] = now

    t0 = time.time()
    stats = demux_to_bins(
        args.r1, panel, out_dir=args.out,
        r2_path=args.r2, policy=policy, workers=workers,
        chunk_size=args.chunk_size,
        output_format=args.format, compress=not args.no_compress,
        split_by_strand=args.split_by_strand,
        progress=progress,
    )
    dt = time.time() - t0
    print(f"\r  done.                              ", file=sys.stderr)

    print(f"\nDemuxed {stats.total:,} reads in {dt:.1f}s "
          f"({stats.total/dt:,.0f} reads/sec)")
    print(f"Assigned: {stats.assigned:,} ({stats.assigned_fraction:.1%})   "
          f"Unassigned: {stats.unassigned:,}")
    print(f"\nPer-bin counts (top 10):")
    for k, v in stats.bin_counts.most_common(10):
        print(f"  {k:20s} {v:>10,}")
    if stats.reason_counts:
        print(f"\nUnassigned reasons:")
        for k, v in stats.reason_counts.most_common():
            print(f"  {k:28s} {v:>10,}")
    print(f"\nOutput written to: {Path(args.out).resolve()}")

    if args.stats_json:
        Path(args.stats_json).write_text(json.dumps(stats.as_dict(), indent=2))
        print(f"Stats JSON written to: {Path(args.stats_json).resolve()}")

    if args.report:
        from ampliscan.report import write_report
        report_path = Path(args.out) / "report.html" if args.report is True else Path(args.report)
        params = {
            "anchor mismatches": policy.anchor_mismatch,
            "barcode mismatches": policy.barcode_mismatch,
            "indels": "yes" if policy.allow_indels else "no",
            "input": Path(args.r1).name + (f" + {Path(args.r2).name}" if args.r2 else ""),
        }
        title = Path(args.r1).name.replace("_R1", "").split(".")[0]
        written = write_report(stats, panel, report_path, title=title, params=params)
        print(f"HTML report written to: {written}")
    return 0


def cmd_validate(args) -> int:
    try:
        panel = _load_panel(args.config)
    except Exception as e:  # noqa: BLE001
        print(f"INVALID: {e}", file=sys.stderr)
        return 1
    print("OK\n")
    print(panel.describe())
    return 0


def cmd_sniff(args) -> int:
    from itertools import islice
    from collections import Counter
    from ampliscan import demux_records
    from ampliscan.io import iter_records

    panel = _load_panel(args.config)
    policy = _make_policy(args)

    r1 = islice(iter_records(args.r1), args.n)
    r2 = islice(iter_records(args.r2), args.n) if args.r2 else None
    results = list(demux_records(r1, panel, policy, mates=r2))
    n = len(results)
    counts = Counter(res.bin_name or f"unassigned:{res.reason}" for res in results)
    assigned = sum(v for k, v in counts.items() if not k.startswith("unassigned"))

    print(f"Previewed {n:,} reads")
    print(f"Assigned: {assigned:,} ({assigned/n:.1%})\n")
    for k, v in counts.most_common(20):
        print(f"  {k:35s} {v:>7,} ({v/n:.1%})")
    return 0


def cmd_gui(args) -> int:
    from ampliscan.gui import launch
    return launch()


def cmd_panel(args) -> int:
    from ampliscan import panel_from_excel

    panel = panel_from_excel(
        args.excel, sheet_name=args.sheet,
        barcode_length=args.barcode_length, name=args.name,
    )
    panel.to_yaml(args.out)
    print(panel.describe())
    print(f"\nPanel written to: {Path(args.out).resolve()}")
    return 0


def _add_match_args(p):
    p.add_argument("--anchor-mismatch", type=int, default=0,
                   help="substitutions allowed in the constant anchor (default 0)")
    p.add_argument("--barcode-mismatch", type=int, default=1,
                   help="substitutions allowed in the barcode (default 1)")
    p.add_argument("--indels", action="store_true",
                   help="use Levenshtein (allow indels) instead of Hamming")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ampliscan",
        description="Anchor-based amplicon demultiplexer.",
    )
    parser.add_argument("--version", action="version",
                        version=f"ampliscan {__version__}")
    sub = parser.add_subparsers(dest="command")

    # gui
    g = sub.add_parser("gui", help="launch the desktop GUI (no terminal needed)")
    g.set_defaults(func=cmd_gui)

    # demux
    d = sub.add_parser("demux", help="demultiplex FASTQ/FASTA into per-bin files")
    d.add_argument("--r1", required=True, help="R1 FASTQ/FASTA (gz ok)")
    d.add_argument("--r2", default=None, help="R2 FASTQ/FASTA for paired-end (gz ok)")
    d.add_argument("--config", required=True, help="panel YAML")
    d.add_argument("--out", required=True, help="output directory for bins")
    d.add_argument("--workers", type=int, default=1,
                   help="CPU cores for matching (default 1; 0=all). Raising this "
                        "only helps for uncompressed input, --indels, or huge panels "
                        "-- for gzipped Illumina input install ampliscan[fast] instead")
    d.add_argument("--chunk-size", type=int, default=20000,
                   help="reads per work chunk (default 20000)")
    d.add_argument("--format", choices=["fastq", "fasta", "both"], default="fastq",
                   help="output format (default fastq)")
    d.add_argument("--no-compress", action="store_true",
                   help="write plain (uncompressed) output")
    d.add_argument("--split-by-strand", action="store_true",
                   help="separate + and - strands into F01_R01_fwd / _rev files")
    d.add_argument("--report", nargs="?", const=True, default=None, metavar="PATH",
                   help="write a self-contained HTML report (heatmaps + unassigned-reason "
                        "chart); PATH defaults to <out>/report.html")
    d.add_argument("--stats-json", default=None,
                   help="also write run stats to this JSON file")
    _add_match_args(d)
    d.set_defaults(func=cmd_demux)

    # validate
    v = sub.add_parser("validate", help="lint a panel YAML")
    v.add_argument("config", help="panel YAML")
    v.set_defaults(func=cmd_validate)

    # sniff
    s = sub.add_parser("sniff", help="preview bin distribution on first N reads")
    s.add_argument("--r1", required=True)
    s.add_argument("--r2", default=None)
    s.add_argument("--config", required=True)
    s.add_argument("-n", type=int, default=10000, help="reads to preview (default 10000)")
    _add_match_args(s)
    s.set_defaults(func=cmd_sniff)

    # panel
    p = sub.add_parser("panel", help="build a panel YAML from an Excel sheet")
    p.add_argument("excel", help="path to the Excel amplicon table")
    p.add_argument("--out", required=True, help="output panel YAML path")
    p.add_argument("--sheet", default=0, help="sheet name or index (default 0)")
    p.add_argument("--barcode-length", type=int, default=4)
    p.add_argument("--name", default="panel_from_excel")
    p.set_defaults(func=cmd_panel)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # No subcommand -> launch the GUI (the friendly default for double-clickers).
    if getattr(args, "command", None) is None:
        return cmd_gui(args)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
