"""Regression tests for (1) the 3' anchor being picked inside the sample index
after the adapter and (2) the failure reason reported for unassigned pairs."""
from ampliscan import Panel, MatchPolicy, demux_record
from ampliscan.anchors import locate_barcoded_region
from ampliscan.io import Record

# The real panel: GCTT / ACAG anchors, 4 nt barcodes.
PANEL = Panel(
    name="real",
    forward_5p_anchor="GCTT",
    forward_3p_anchor="ACAG",
    forward_barcodes={"F01": "GCGT", "F02": "GTAG", "F03": "ACGC", "F04": "CTCG", "F05": "GCTC"},
    reverse_barcodes={"R01": "ACGC", "R02": "CTAC", "R03": "GCGT", "R04": "CGAG"},
)

TARGET = "ATGGATACGCATTCAAAGCCTTGGAGCCATTG"
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"      # Illumina adapter up to the i7 index
INDEX_WITH_ACAG = "AGTATCTTGC"                       # "...CCAGTC ACAG TATC..." has a 2nd ACAG
INDEX_NO_ACAG = "TTACGTTCGA"
TAIL = "ATCTCGTATGCCG"


def _read(index, target=TARGET, f="GCGT", r="ACGC", adapter=ADAPTER):
    return "GCTT" + f + target + r + "ACAG" + adapter + index + TAIL


def _rec(seq, rid="r"):
    return Record(id=rid, seq=seq, qual="I" * len(seq))


def test_stray_anchor_in_index_is_not_used():
    # Exact read from the bug report. Without trimming, the rightmost ACAG is
    # the one inside "CCAGTC ACAG TATC" (the index), not the real anchor.
    seq = ("GCTTGCGT" + TARGET + "ACGC" + "ACAG"
           + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTATCTTGCATCTCGTATGCCG")
    res = demux_record(_rec(seq), PANEL)
    assert res.bin_name == "F01_R01"
    assert res.target_seq == TARGET


def test_index_without_acag_still_assigns_identically():
    a = demux_record(_rec(_read(INDEX_WITH_ACAG)), PANEL)
    b = demux_record(_rec(_read(INDEX_NO_ACAG)), PANEL)
    assert a.bin_name == b.bin_name == "F01_R01"
    assert a.target_seq == b.target_seq == TARGET


def test_read_without_adapter_is_unchanged_and_rightmost_anchor_wins():
    # No adapter, and an ACAG inside the target: the rightmost ACAG (the real one)
    # must still win, and the result must equal the untrimmed search.
    target = "AAACAGGGTTTGACTATCACC"
    seq = "GCTT" + "GCGT" + target + "CGAG" + "ACAG" + "TTTTTTTT"
    on = demux_record(_rec(seq), PANEL)
    off = demux_record(_rec(seq), PANEL, MatchPolicy(adapter_trim=""))
    assert on.bin_name == off.bin_name == "F01_R04"
    assert on.target_seq == off.target_seq == target


def test_internal_anchor_before_adapter_still_rightmost_within_window():
    target = "AAACAGGGTTTGACTATCACC"
    seq = _read(INDEX_WITH_ACAG, target=target, r="CTAC")
    res = demux_record(_rec(seq), PANEL)
    assert res.bin_name == "F01_R02"
    assert res.target_seq == target


def test_trimming_disabled_gives_old_behaviour():
    seq = _read(INDEX_WITH_ACAG)
    res = demux_record(_rec(seq), PANEL, MatchPolicy(adapter_trim=""))
    assert res.bin_name is None and res.reason == "barcode_unmatched"
    # ... and the low-level function agrees with the policy switch
    old = locate_barcoded_region(seq, "GCTT", "ACAG", 4)
    new = locate_barcoded_region(seq, "GCTT", "ACAG", 4, adapter="AGATCGGAAGAGC")
    assert old[3] != TARGET and new[3] == TARGET


def test_target_never_contains_adapter_bases():
    res = demux_record(_rec(_read(INDEX_WITH_ACAG)), PANEL)
    assert "AGATCGGAAGAGC" not in res.target_seq


def test_adapter_before_5p_barcode_is_ignored():
    # Adapter text earlier than the F barcode must not shrink the window.
    seq = "AGATCGGAAGAGC" + _read(INDEX_NO_ACAG)
    assert demux_record(_rec(seq), PANEL).bin_name == "F01_R01"


def test_reverse_strand_read_is_trimmed_too():
    from ampliscan.config import reverse_complement
    # An R2-like read is the reverse complement of the molecule; the adapter is
    # still at its 3' end. Build the molecule with the index-borne ACAG, then
    # make the read as the sequencer would: RC molecule + adapter + index.
    mol = "GCTT" + "GCGT" + TARGET + "ACGC" + "ACAG"
    rc = reverse_complement(mol)
    # The index text contains a stray AAGC (the reverse-strand 3' anchor).
    seq = rc + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGT" + "TTAAGCTTACGT"
    res = demux_record(_rec(seq), PANEL)
    assert res.bin_name == "F01_R01" and res.strand == "-"
    assert res.target_seq == TARGET
    # and without trimming the stray AAGC is taken instead
    assert demux_record(_rec(seq), PANEL, MatchPolicy(adapter_trim="")).bin_name is None


# ---- Bug 2: reason reported for unassigned pairs -------------------------------

def _n_after_8(n=151):
    return "ACGTACGT" + "N" * (n - 8)


def test_pair_unassigned_reports_informative_reason_not_no_anchor():
    r1 = _read(INDEX_NO_ACAG, r="TTTT")          # anchors present, R barcode unknown
    res = demux_record(_rec(r1), PANEL, mate=_rec(_n_after_8()))
    assert res.bin_name is None
    assert res.reason == "barcode_unmatched"


def test_pair_assigned_by_r1_when_r2_is_all_n():
    res = demux_record(_rec(_read(INDEX_NO_ACAG)), PANEL, mate=_rec(_n_after_8()))
    assert res.bin_name == "F01_R01"
    assert res.r1_bin_name == "F01_R01" and res.r2_bin_name is None


def test_pair_assigned_by_r2_when_r1_has_no_anchor():
    from ampliscan.config import reverse_complement
    mol = "GCTT" + "GCGT" + TARGET + "ACGC" + "ACAG"
    r2 = reverse_complement(mol) + "AGATCGGAAGAGCGTCGTGTAGGG"
    res = demux_record(_rec(_n_after_8()), PANEL, mate=_rec(r2))
    assert res.bin_name == "F01_R01"
    assert res.r1_bin_name is None and res.r2_bin_name == "F01_R01"


def test_pair_reason_comes_from_r2_when_r1_has_no_anchor():
    from ampliscan.config import reverse_complement
    bad = "GCTT" + "GCGT" + TARGET + "TTTT" + "ACAG"
    r2 = reverse_complement(bad) + "AGATCGGAAGAGCGTCGTGTAGGG"
    res = demux_record(_rec(_n_after_8()), PANEL, mate=_rec(r2))
    assert res.bin_name is None and res.reason == "barcode_unmatched"


def test_pair_both_no_anchor_reports_no_anchor():
    res = demux_record(_rec(_n_after_8()), PANEL, mate=_rec(_n_after_8()))
    assert res.bin_name is None and res.reason == "no_anchor"


def test_pair_disagreement_unchanged():
    r1 = _read(INDEX_NO_ACAG, f="GCGT", r="ACGC")   # F01_R01
    r2 = _read(INDEX_NO_ACAG, f="GTAG", r="ACGC")   # F02_R01
    res = demux_record(_rec(r1), PANEL, mate=_rec(r2))
    assert res.reason == "paired_end_disagreement"
