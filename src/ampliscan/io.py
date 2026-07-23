"""Sequence record I/O: FASTQ + FASTA, gzip-aware."""
from __future__ import annotations

import gzip
import io as _io
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

# Prefer isal (igzip) when available: 2-5x faster gzip read/write, drop-in.
# The real bottleneck for gzipped Illumina input is (de)compression, not the
# barcode matching, so this is the single highest-value speedup.
try:  # pragma: no cover - optional accelerator
    from isal import igzip as _fastgzip
    HAVE_ISAL = True
except ImportError:  # pragma: no cover
    _fastgzip = gzip
    HAVE_ISAL = False


@dataclass
class Record:
    id: str
    seq: str
    qual: Optional[str] = None      # None for FASTA
    description: str = ""

    @property
    def is_fastq(self) -> bool:
        return self.qual is not None


PathLike = Union[str, Path]


def _open(path: PathLike, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        # isal's igzip.open writes a standard gzip stream readable by any tool.
        # For writing, compresslevel is capped at 3 in isal but is much faster
        # and produces comparable sizes for FASTQ text.
        return _fastgzip.open(path, mode)
    return open(path, mode)


def _detect_format(path: PathLike) -> str:
    name = str(path).lower().rstrip(".gz")
    if name.endswith((".fq", ".fastq")):
        return "fastq"
    if name.endswith((".fa", ".fasta", ".fna")):
        return "fasta"
    raise ValueError(f"Cannot infer format from extension: {path}")


def iter_fastq(path: PathLike) -> Iterator[Record]:
    with _open(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                raise ValueError(f"Truncated FASTQ record near: {header!r}")
            header = header.rstrip("\n")
            if not header.startswith("@"):
                raise ValueError(f"Bad FASTQ header (no @): {header!r}")
            header = header[1:]
            rid, _, desc = header.partition(" ")
            yield Record(id=rid, seq=seq, qual=qual, description=desc)


def iter_fasta(path: PathLike) -> Iterator[Record]:
    with _open(path, "rt") as fh:
        rid = None
        desc = ""
        chunks: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rid is not None:
                    yield Record(id=rid, seq="".join(chunks), qual=None, description=desc)
                header = line[1:]
                rid, _, desc = header.partition(" ")
                chunks = []
            else:
                chunks.append(line)
        if rid is not None:
            yield Record(id=rid, seq="".join(chunks), qual=None, description=desc)


def iter_records(path: PathLike, format: Optional[str] = None) -> Iterator[Record]:
    """Iterate records from FASTQ or FASTA (auto-detected by extension)."""
    fmt = format or _detect_format(path)
    if fmt == "fastq":
        yield from iter_fastq(path)
    elif fmt == "fasta":
        yield from iter_fasta(path)
    else:
        raise ValueError(f"Unknown format: {fmt}")


def write_fastq(records: Iterable[Record], path: PathLike) -> int:
    n = 0
    with _open(path, "wt") as fh:
        for r in records:
            qual = r.qual if r.qual is not None else "I" * len(r.seq)
            desc = f" {r.description}" if r.description else ""
            fh.write(f"@{r.id}{desc}\n{r.seq}\n+\n{qual}\n")
            n += 1
    return n


def write_fasta(records: Iterable[Record], path: PathLike, line_width: int = 80) -> int:
    n = 0
    with _open(path, "wt") as fh:
        for r in records:
            desc = f" {r.description}" if r.description else ""
            fh.write(f">{r.id}{desc}\n")
            seq = r.seq
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + "\n")
            n += 1
    return n


def write_bins(
    results: Iterable,
    out_dir: PathLike,
    output_format: str = "fastq",
    compress: bool = False,
    include_unassigned: bool = True,
) -> dict[str, int]:
    """Group DemuxResults by bin_name and write one file per bin.

    Returns a dict mapping bin_name -> read count.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ext = "fq" if output_format == "fastq" else "fa"
    if compress:
        ext += ".gz"

    buckets: dict[str, list[Record]] = defaultdict(list)
    for res in results:
        bin_name = res.bin_name or "_unassigned"
        if bin_name == "_unassigned" and not include_unassigned:
            continue
        seq = res.target_seq if res.target_seq is not None else ""
        qual = res.target_qual if res.target_qual is not None else (
            "I" * len(seq) if output_format == "fastq" else None
        )
        rec = Record(id=res.read_id, seq=seq, qual=qual,
                     description=f"strand={res.strand} reason={res.reason or 'ok'}")
        buckets[bin_name].append(rec)

    counts = {}
    for bin_name, recs in buckets.items():
        path = out_dir / f"{bin_name}.{ext}"
        if output_format == "fastq":
            counts[bin_name] = write_fastq(recs, path)
        else:
            counts[bin_name] = write_fasta(recs, path)
    return counts
