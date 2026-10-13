"""Regression: reads with no insert between the barcodes must never reach a sample bin.

Real S25 data contained primer-dimer-like reads such as
``GCTT GCGT GAGT ACAG AGATCGGAAGAG...`` (F barcode runs straight into the 3' anchor;
``GAGT`` is 1 mismatch from R03 ``GCGT``). They were written as empty FASTQ
records, and one empty record crashes CRISPResso2's aligner for the whole sample.
"""
import gzip

from ampliscan import MatchPolicy, Panel, demux_record, demux_to_bins
from ampliscan.cli import main
from ampliscan.io import Record, write_fastq

PANEL = Panel(
    name="t",
    forward_5p_anchor="GCTT",
    forward_3p_anchor="ACAG",
    forward_barcodes={"F01": "GCGT", "F02": "GTAG"},
    reverse_barcodes={"R01": "ACGC", "R03": "GCGT"},
)
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
DIMER = "GCTT" + "GCGT" + "GAGT" + "ACAG" + ADAPTER          # insert is empty
GOOD = "GCTT" + "GCGT" + "TTGACCATGGCCAATTGGCA" + "GCGT" + "ACAG" + ADAPTER


def test_empty_insert_is_not_assigned():
    res = demux_record(Record("d", DIMER, "I" * len(DIMER)), PANEL, MatchPolicy(barcode_mismatch=1))
    assert res.bin_name is None
    assert res.reason == "target_too_short"
    assert res.target_seq == ""


def test_normal_read_still_assigned():
    res = demux_record(Record("g", GOOD, "I" * len(GOOD)), PANEL, MatchPolicy(barcode_mismatch=1))
    assert res.bin_name == "F01_R03"
    assert len(res.target_seq) == 20


def test_min_target_length_is_configurable():
    rec = Record("g", GOOD, "I" * len(GOOD))
    assert demux_record(rec, PANEL, MatchPolicy(min_target_length=21)).reason == "target_too_short"
    assert demux_record(rec, PANEL, MatchPolicy(min_target_length=20)).bin_name == "F01_R03"


def test_bins_contain_no_empty_records(tmp_path):
    recs = [Record(f"g{i}", GOOD, "I" * len(GOOD)) for i in range(30)]
    recs += [Record(f"d{i}", DIMER, "I" * len(DIMER)) for i in range(5)]
    p = tmp_path / "r1.fq.gz"
    write_fastq(recs, p)

    stats = demux_to_bins(p, PANEL, tmp_path / "out", policy=MatchPolicy(barcode_mismatch=1))
    assert stats.assigned == 30
    assert stats.reason_counts["target_too_short"] == 5

    with gzip.open(tmp_path / "out" / "F01_R03.fq.gz", "rt") as fh:
        seqs = [l.rstrip("\n") for i, l in enumerate(fh) if i % 4 == 1]
    assert len(seqs) == 30 and all(seqs)             # nothing empty in the bin


def test_cli_min_target_length_flag(tmp_path):
    recs = [Record(f"g{i}", GOOD, "I" * len(GOOD)) for i in range(10)]
    p = tmp_path / "r1.fq.gz"
    write_fastq(recs, p)
    cfg = tmp_path / "panel.yaml"
    PANEL.to_yaml(cfg)
    rc = main(["demux", "--r1", str(p), "--config", str(cfg), "--out", str(tmp_path / "o"),
               "--min-target-length", "50"])
    assert rc == 0
    assert not (tmp_path / "o" / "F01_R03.fq.gz").exists()   # all 20-nt inserts fall below 50
