"""Smoke tests covering the end-to-end synthetic pipeline."""
from ampliscan import (
    Panel, MatchPolicy, simulate_reads, demux_records, panel_from_excel,
)
from ampliscan.matcher import BarcodeMatcher, hamming, levenshtein
from ampliscan.anchors import find_anchor, locate_barcoded_region


SMALL_PANEL = Panel(
    name="test",
    forward_5p_anchor="GCTT",
    forward_3p_anchor="ACAG",
    forward_barcodes={"F01": "GCGT", "F02": "GTAG", "F03": "ACGC"},
    reverse_barcodes={"R01": "ACGC", "R02": "CTAC"},
)


def test_panel_inference():
    assert SMALL_PANEL.barcode_length == 4
    assert SMALL_PANEL.reverse_5p_anchor == "CTGT"  # revcomp(ACAG)
    assert SMALL_PANEL.reverse_3p_anchor == "AAGC"  # revcomp(GCTT)


def test_hamming_and_levenshtein():
    assert hamming("ACGT", "ACGT") == 0
    assert hamming("ACGT", "ACCT") == 1
    assert hamming("ACGT", "TCCA") == 3
    assert levenshtein("ACGT", "AGT") == 1
    assert levenshtein("ACGT", "ACGTACGT", max_dist=2) == 3  # exceeds budget


def test_barcode_matcher_exact_and_fuzzy():
    m = BarcodeMatcher({"A": "GCGT", "B": "GTAG"}, max_mismatch=1)
    assert m.match("GCGT") == ("A", 0, False)
    assert m.match("GCAT") == ("A", 1, False)
    assert m.match("GCAA")[0] is None  # too far


def test_find_anchor_picks_rightmost_with_from_end():
    read = "AAAGCTTNNNNGCTTTTTT"
    assert find_anchor(read, "GCTT") == 3
    assert find_anchor(read, "GCTT", from_end=True) == 11


def test_locate_barcoded_region():
    seq = "XXGCTT" + "GCGT" + "TARGETSEQ" + "ACGC" + "ACAG" + "YY"
    hit = locate_barcoded_region(seq, "GCTT", "ACAG", barcode_length=4)
    assert hit is not None
    _, _, fbc, target, rbc = hit
    assert fbc == "GCGT"
    assert rbc == "ACGC"
    assert target == "TARGETSEQ"


def test_end_to_end_perfect_data():
    r1, r2 = simulate_reads(SMALL_PANEL, reads_per_bin=20, target_length=80,
                            error_rate=0.0, paired_end=True, seed=0)
    policy = MatchPolicy(anchor_mismatch=0, barcode_mismatch=0)
    results = list(demux_records(iter(r1), SMALL_PANEL, policy, mates=iter(r2)))
    correct = 0
    for res in results:
        parts = res.read_id.split("_")
        true_bin = f"{parts[1]}_{parts[2]}"
        if res.bin_name == true_bin:
            correct += 1
    assert correct == len(results), f"only {correct}/{len(results)} correct"


def test_end_to_end_noisy_data_high_accuracy():
    r1, r2 = simulate_reads(SMALL_PANEL, reads_per_bin=100, target_length=100,
                            error_rate=0.005, paired_end=True, seed=1)
    policy = MatchPolicy(anchor_mismatch=0, barcode_mismatch=1)
    results = list(demux_records(iter(r1), SMALL_PANEL, policy, mates=iter(r2)))
    assigned = sum(1 for r in results if r.bin_name)
    assert assigned / len(results) >= 0.95


def test_r1_r2_bin_names_populated_and_balanced():
    """Each mate's own bin assignment must be tracked independently of the
    paired reconciliation, and (for realistic paired data) roughly balanced --
    not skewed toward R1 just because R1 is checked first."""
    r1, r2 = simulate_reads(SMALL_PANEL, reads_per_bin=150, target_length=90,
                            error_rate=0.005, paired_end=True, seed=3)
    policy = MatchPolicy(anchor_mismatch=0, barcode_mismatch=1)
    results = list(demux_records(iter(r1), SMALL_PANEL, policy, mates=iter(r2)))

    n_r1 = sum(1 for res in results if res.r1_bin_name)
    n_r2 = sum(1 for res in results if res.r2_bin_name)
    assert n_r1 > 0 and n_r2 > 0
    # neither mate should dominate -- within 10% of each other
    assert abs(n_r1 - n_r2) / max(n_r1, n_r2) < 0.10, (n_r1, n_r2)


def test_r1_r2_bin_names_agree_when_pair_assigned():
    r1, r2 = simulate_reads(SMALL_PANEL, reads_per_bin=50, target_length=90,
                            error_rate=0.0, paired_end=True, seed=4)
    results = list(demux_records(iter(r1), SMALL_PANEL,
                                 MatchPolicy(barcode_mismatch=0), mates=iter(r2)))
    for res in results:
        if res.bin_name:  # pair was assigned -> both mates must agree with it
            assert res.r1_bin_name == res.bin_name
            assert res.r2_bin_name == res.bin_name


def test_r1_r2_bin_names_single_end():
    """Single-end input has no mate, so r2_bin_name must stay None."""
    r1, _ = simulate_reads(SMALL_PANEL, reads_per_bin=20, target_length=80,
                           error_rate=0.0, paired_end=False, seed=5)
    results = list(demux_records(iter(r1), SMALL_PANEL, MatchPolicy(barcode_mismatch=0)))
    assert all(res.r2_bin_name is None for res in results)
    assert any(res.r1_bin_name is not None for res in results)


def test_excel_panel_parses_real_file():
    import os
    import pytest
    xlsx = r"C:/Users/ayush/Downloads/Data_Ayushman.xlsx"
    if not os.path.exists(xlsx):
        pytest.skip("Data_Ayushman.xlsx not present on this machine")
    p = panel_from_excel(xlsx, barcode_length=4)
    assert p.forward_5p_anchor == "GCTT"
    assert p.forward_3p_anchor == "ACAG"
    assert p.forward_barcodes["F01"] == "GCGT"
    assert p.reverse_barcodes["R01"] == "ACGC"
