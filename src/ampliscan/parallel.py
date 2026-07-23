"""Parallel demultiplexing across multiple CPU cores.

Reads FASTQ (optionally paired) in chunks, farms each chunk out to a worker
process, and streams the returned results to per-bin output files so peak
memory stays bounded even on multi-million-read sequencer files.

Windows-safe: workers use the ``spawn`` start method, so the worker function
and its state initializer are module-level and picklable. Callers that invoke
:func:`demux_to_bins` from a script must guard the entry point with
``if __name__ == "__main__":``.
"""
from __future__ import annotations

import gzip
import os
from collections import Counter
from dataclasses import dataclass, field
from itertools import islice
from multiprocessing import get_context
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple, Union

from ampliscan.config import Panel, MatchPolicy
from ampliscan.core import _demux_with_matchers, DemuxResult
from ampliscan.matcher import BarcodeMatcher
from ampliscan.io import Record, iter_records, PathLike, _fastgzip


# --- worker-side state (one per process, built once via initializer) ---------

_WORKER: dict = {}


def _init_worker(panel: Panel, policy: MatchPolicy) -> None:
    _WORKER["panel"] = panel
    _WORKER["policy"] = policy
    _WORKER["fwd"] = BarcodeMatcher(
        panel.forward_barcodes, policy.barcode_mismatch, policy.allow_indels
    )
    _WORKER["rev"] = BarcodeMatcher(
        panel.reverse_barcodes, policy.barcode_mismatch, policy.allow_indels
    )


def _process_chunk(
    chunk: List[Tuple[Record, Optional[Record]]]
) -> List[DemuxResult]:
    panel = _WORKER["panel"]
    policy = _WORKER["policy"]
    fwd = _WORKER["fwd"]
    rev = _WORKER["rev"]
    out = []
    for r1, r2 in chunk:
        out.append(_demux_with_matchers(r1, panel, policy, fwd, rev, r2))
    return out


# --- chunking helpers --------------------------------------------------------

def _paired_iter(
    r1_path: PathLike, r2_path: Optional[PathLike], fmt: Optional[str]
) -> Iterator[Tuple[Record, Optional[Record]]]:
    it1 = iter_records(r1_path, format=fmt)
    if r2_path is None:
        for rec in it1:
            yield rec, None
    else:
        it2 = iter_records(r2_path, format=fmt)
        for a, b in zip(it1, it2):
            yield a, b


def _chunked(
    it: Iterator[Tuple[Record, Optional[Record]]], size: int
) -> Iterator[List[Tuple[Record, Optional[Record]]]]:
    while True:
        block = list(islice(it, size))
        if not block:
            return
        yield block


# --- output stats ------------------------------------------------------------

@dataclass
class DemuxStats:
    total: int = 0
    assigned: int = 0
    bin_counts: Counter = field(default_factory=Counter)
    reason_counts: Counter = field(default_factory=Counter)
    # Per-mate breakdown: what bin would R1 alone assign, and what bin would
    # R2 alone assign, independent of the paired-end reconciliation. Empty
    # for single-end input. Useful as a QC heatmap comparing whether R1 and
    # R2 independently agree on per-sample counts.
    bin_counts_r1: Counter = field(default_factory=Counter)
    bin_counts_r2: Counter = field(default_factory=Counter)
    workers: int = 1

    @property
    def unassigned(self) -> int:
        return self.total - self.assigned

    @property
    def assigned_fraction(self) -> float:
        return self.assigned / self.total if self.total else 0.0

    def as_dict(self) -> dict:
        return {
            "total": self.total,
            "assigned": self.assigned,
            "unassigned": self.unassigned,
            "assigned_fraction": self.assigned_fraction,
            "workers": self.workers,
            "bin_counts": dict(self.bin_counts),
            "bin_counts_r1": dict(self.bin_counts_r1),
            "bin_counts_r2": dict(self.bin_counts_r2),
            "reason_counts": dict(self.reason_counts),
        }


# --- streaming per-bin writer ------------------------------------------------

_EXT = {"fastq": "fq", "fasta": "fa"}


def _normalize_formats(output_format) -> list:
    """Accept 'fastq' | 'fasta' | 'both' | a list, return an ordered list."""
    if isinstance(output_format, (list, tuple, set)):
        fmts = list(output_format)
    elif output_format == "both":
        fmts = ["fastq", "fasta"]
    else:
        fmts = [output_format]
    for f in fmts:
        if f not in _EXT:
            raise ValueError(f"Unknown output format: {f!r} (use fastq/fasta/both)")
    # de-dupe, keep order
    seen, out = set(), []
    for f in fmts:
        if f not in seen:
            seen.add(f); out.append(f)
    return out


class _BinWriter:
    """Lazily-opened per-bin output files, written as results stream in.

    Can emit more than one format at once (e.g. FASTQ *and* FASTA) -- each
    (bin, format) pair gets its own file, opened on first use.
    """

    def __init__(self, out_dir: Path, output_format, compress: bool,
                 split_by_strand: bool = False):
        self.out_dir = out_dir
        self.formats = _normalize_formats(output_format)
        self.compress = compress
        self.split_by_strand = split_by_strand
        self._handles: Dict[tuple, object] = {}
        out_dir.mkdir(parents=True, exist_ok=True)

    def _handle(self, bin_name: str, fmt: str):
        key = (bin_name, fmt)
        h = self._handles.get(key)
        if h is None:
            ext = _EXT[fmt] + (".gz" if self.compress else "")
            path = self.out_dir / f"{bin_name}.{ext}"
            h = _fastgzip.open(path, "wt") if self.compress else open(path, "wt")
            self._handles[key] = h
        return h

    def write(self, res: DemuxResult) -> None:
        bin_name = res.bin_name or "_unassigned"
        # Optionally separate the + strand (matched forward anchors, R1-like)
        # from the - strand (matched revcomp anchors, R2-like) into own files.
        if self.split_by_strand and res.bin_name:
            suffix = {"+": "fwd", "-": "rev"}.get(res.strand, "na")
            bin_name = f"{bin_name}_{suffix}"
        seq = res.target_seq if res.target_seq is not None else ""
        desc = f"strand={res.strand} reason={res.reason or 'ok'}"
        for fmt in self.formats:
            h = self._handle(bin_name, fmt)
            if fmt == "fastq":
                qual = res.target_qual if res.target_qual is not None else "I" * len(seq)
                h.write(f"@{res.read_id} {desc}\n{seq}\n+\n{qual}\n")
            else:
                h.write(f">{res.read_id} {desc}\n{seq}\n")

    def close(self) -> None:
        for h in self._handles.values():
            h.close()
        self._handles.clear()


# --- public entry point ------------------------------------------------------

def demux_to_bins(
    r1_path: PathLike,
    panel: Panel,
    out_dir: PathLike,
    r2_path: Optional[PathLike] = None,
    policy: Optional[MatchPolicy] = None,
    workers: int = 1,
    chunk_size: int = 20000,
    output_format: str = "fastq",
    compress: bool = True,
    split_by_strand: bool = False,
    format: Optional[str] = None,
    progress: Optional[callable] = None,
) -> DemuxStats:
    """Demultiplex FASTQ/FASTA files to per-bin output.

    ``workers`` controls how many CPU cores run the barcode matching in
    parallel. **The default is 1 (single process) on purpose.**

    For typical gzipped Illumina input the job is I/O-bound (gzip
    decompression/compression dominates, not the matching), so spreading the
    matching across cores does not help and the inter-process pickling of
    millions of records actually makes it *slower*. The real speedup for that
    case is installing ``isal`` (``pip install ampliscan[fast]``), which
    accelerates gzip 2-5x and is used automatically when present.

    Raise ``workers`` above 1 only when the matching itself is the bottleneck:
    uncompressed input, Levenshtein matching (``allow_indels=True``), or very
    large barcode panels. Pass ``workers=0`` to auto-detect all cores.

    ``progress`` if given is called as ``progress(n_done)`` after each chunk.
    Returns a :class:`DemuxStats` summarising per-bin counts and rejection
    reasons. Output is written incrementally, so peak memory is roughly
    ``chunk_size * workers`` reads regardless of total input size.
    """
    policy = policy or MatchPolicy()
    if workers == 0:
        workers = os.cpu_count() or 1
    workers = max(1, int(workers))

    out_dir = Path(out_dir)
    writer = _BinWriter(out_dir, output_format, compress, split_by_strand=split_by_strand)
    stats = DemuxStats(workers=workers)

    def _tally_and_write(results: Iterable[DemuxResult]) -> None:
        for res in results:
            stats.total += 1
            if res.bin_name:
                stats.assigned += 1
                stats.bin_counts[res.bin_name] += 1
            else:
                stats.reason_counts[res.reason or "unknown"] += 1
            # Per-mate tallies are independent of the overall pair outcome --
            # R1 (or R2) can assign on its own even when the pair disagrees
            # or the other mate fails, and that's exactly what the QC
            # heatmap should surface.
            if res.r1_bin_name:
                stats.bin_counts_r1[res.r1_bin_name] += 1
            if res.r2_bin_name:
                stats.bin_counts_r2[res.r2_bin_name] += 1
            writer.write(res)

    try:
        pairs = _paired_iter(r1_path, r2_path, format)
        chunks = _chunked(pairs, chunk_size)

        if workers == 1:
            _init_worker(panel, policy)
            for chunk in chunks:
                _tally_and_write(_process_chunk(chunk))
                if progress:
                    progress(stats.total)
        else:
            ctx = get_context("spawn")
            with ctx.Pool(
                processes=workers,
                initializer=_init_worker,
                initargs=(panel, policy),
            ) as pool:
                # imap keeps output order and bounds in-flight chunks
                for results in pool.imap(_process_chunk, chunks):
                    _tally_and_write(results)
                    if progress:
                        progress(stats.total)
    finally:
        writer.close()

    return stats
