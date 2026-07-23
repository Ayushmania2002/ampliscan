"""Tests for the file-based / parallel demux path (demux_to_bins)."""
import gzip
from pathlib import Path

import pytest

from ampliscan import Panel, MatchPolicy, simulate_reads, demux_to_bins
from ampliscan.io import write_fastq


PANEL = Panel(
    name="test",
    forward_5p_anchor="GCTT",
    forward_3p_anchor="ACAG",
    forward_barcodes={"F01": "GCGT", "F02": "GTAG", "F03": "ACGC"},
    reverse_barcodes={"R01": "ACGC", "R02": "CTAC"},
)


@pytest.fixture
def paired_files(tmp_path):
    r1, r2 = simulate_reads(PANEL, reads_per_bin=50, target_length=80,
                            error_rate=0.0, paired_end=True, seed=0)
    p1 = tmp_path / "r1.fq.gz"
    p2 = tmp_path / "r2.fq.gz"
    write_fastq(r1, p1)
    write_fastq(r2, p2)
    return p1, p2


def test_demux_to_bins_single_process(paired_files, tmp_path):
    p1, p2 = paired_files
    out = tmp_path / "out"
    stats = demux_to_bins(p1, PANEL, out, r2_path=p2,
                          policy=MatchPolicy(barcode_mismatch=0),
                          workers=1, chunk_size=37, compress=True)
    # 3 F x 2 R x 50 reads = 300 pairs, all perfect -> all assigned
    assert stats.total == 300
    assert stats.assigned == 300
    assert stats.assigned_fraction == 1.0
    # one file per bin should exist
    for f in PANEL.forward_barcodes:
        for r in PANEL.reverse_barcodes:
            assert (out / f"{f}_{r}.fq.gz").exists()


def test_demux_to_bins_split_by_strand(paired_files, tmp_path):
    """With split_by_strand, each bin becomes _fwd (+) and _rev (-) files, and
    each file holds only reads of that strand."""
    import gzip
    p1, p2 = paired_files
    out = tmp_path / "out"
    stats = demux_to_bins(p1, PANEL, out, r2_path=p2,
                          policy=MatchPolicy(barcode_mismatch=0),
                          workers=1, split_by_strand=True, compress=True)
    assert stats.total == 300
    # plain (unsplit) bin files must NOT exist
    assert not (out / "F01_R01.fq.gz").exists()
    # at least one strand file exists and is strand-pure
    seen = False
    for f in PANEL.forward_barcodes:
        for r in PANEL.reverse_barcodes:
            for strand, sym in (("fwd", "+"), ("rev", "-")):
                path = out / f"{f}_{r}_{strand}.fq.gz"
                if not path.exists():
                    continue
                seen = True
                with gzip.open(path, "rt") as fh:
                    for i, line in enumerate(fh):
                        if i % 4 == 0 and "strand=" in line:
                            assert f"strand={sym}" in line
    assert seen


def test_demux_to_bins_both_formats(paired_files, tmp_path):
    """output_format='both' writes a .fq.gz and a .fa.gz per bin."""
    p1, p2 = paired_files
    out = tmp_path / "out"
    demux_to_bins(p1, PANEL, out, r2_path=p2,
                  policy=MatchPolicy(barcode_mismatch=0),
                  workers=1, output_format="both", compress=True)
    assert (out / "F01_R01.fq.gz").exists()
    assert (out / "F01_R01.fa.gz").exists()
    # FASTA has 2 lines/record, FASTQ has 4 -> same record count
    import gzip
    with gzip.open(out / "F01_R01.fq.gz", "rt") as fh:
        fq = sum(1 for _ in fh) // 4
    with gzip.open(out / "F01_R01.fa.gz", "rt") as fh:
        fa = sum(1 for _ in fh) // 2
    assert fq == fa and fq > 0


def test_normalize_formats():
    from ampliscan.parallel import _normalize_formats
    assert _normalize_formats("fastq") == ["fastq"]
    assert _normalize_formats("both") == ["fastq", "fasta"]
    assert _normalize_formats(["fasta", "fastq"]) == ["fasta", "fastq"]
    assert _normalize_formats(["fastq", "fastq"]) == ["fastq"]
    import pytest as _pytest
    with _pytest.raises(ValueError):
        _normalize_formats("bam")


def test_demux_to_bins_multiprocess_matches_single(paired_files, tmp_path):
    p1, p2 = paired_files
    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    policy = MatchPolicy(barcode_mismatch=1)
    s1 = demux_to_bins(p1, PANEL, out1, r2_path=p2, policy=policy,
                       workers=1, chunk_size=40, compress=True)
    s2 = demux_to_bins(p1, PANEL, out2, r2_path=p2, policy=policy,
                       workers=2, chunk_size=40, compress=True)
    # Parallelism must not change the result, only how it's computed.
    assert s1.total == s2.total
    assert s1.assigned == s2.assigned
    assert s1.bin_counts == s2.bin_counts
    assert s1.reason_counts == s2.reason_counts


def test_demux_to_bins_chunk_boundary(paired_files, tmp_path):
    """A chunk_size that doesn't divide the total must not drop/dup reads."""
    p1, p2 = paired_files
    out = tmp_path / "out"
    stats = demux_to_bins(p1, PANEL, out, r2_path=p2, workers=1, chunk_size=7)
    assert stats.total == 300


def test_output_files_are_valid_gzip_fastq(paired_files, tmp_path):
    p1, p2 = paired_files
    out = tmp_path / "out"
    demux_to_bins(p1, PANEL, out, r2_path=p2,
                  policy=MatchPolicy(barcode_mismatch=0), workers=1)
    f = out / "F01_R01.fq.gz"
    with gzip.open(f, "rt") as fh:
        lines = fh.readlines()
    assert len(lines) % 4 == 0          # well-formed FASTQ
    assert lines[0].startswith("@")     # header
    assert lines[2].startswith("+")     # separator
