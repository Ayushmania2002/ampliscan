"""Tests for the HTML report and the `demux --report` flag."""
from collections import Counter

from ampliscan import Panel, MatchPolicy, simulate_reads, demux_to_bins
from ampliscan.cli import main
from ampliscan.io import write_fastq
from ampliscan.parallel import DemuxStats
from ampliscan.report import build_report_html, write_report

PANEL = Panel(
    name="test",
    forward_5p_anchor="GCTT",
    forward_3p_anchor="ACAG",
    forward_barcodes={"F01": "GCGT", "F02": "GTAG"},
    reverse_barcodes={"R01": "ACGC", "R02": "CTAC"},
)


def _stats(paired=True, unassigned=True):
    s = DemuxStats(total=1000, assigned=900 if unassigned else 1000)
    s.bin_counts = Counter({"F01_R01": 400, "F01_R02": 200, "F02_R01": 200, "F02_R02": 100})
    if not unassigned:
        s.bin_counts["F02_R02"] = 200
        s.bin_counts["F02_R01"] = 200
        s.bin_counts["F01_R01"] = 400
        s.bin_counts["F01_R02"] = 200
    if unassigned:
        s.reason_counts = Counter({"no_anchor": 60, "barcode_unmatched": 40})
    if paired:
        s.bin_counts_r1 = Counter({"F01_R01": 390, "F02_R02": 110})
        s.bin_counts_r2 = Counter({"F01_R01": 395, "F02_R02": 105})
    return s


def test_report_contains_core_sections():
    page = build_report_html(_stats(), PANEL, title="demo")
    assert "<svg" in page                      # donut chart
    assert 'class="heat"' in page              # heatmap
    assert "No anchor found" in page
    assert "R1 vs R2" in page                  # paired-end per-mate section
    assert "http://" not in page and "https://" not in page.replace("doi:", "")  # no external requests


def test_report_single_end_omits_mate_section():
    page = build_report_html(_stats(paired=False), PANEL)
    assert "R1 vs R2" not in page


def test_report_when_everything_assigned():
    page = build_report_html(_stats(unassigned=False), PANEL)
    assert "Every read was assigned" in page


def test_report_escapes_html_in_title():
    page = build_report_html(_stats(), PANEL, title="<script>alert(1)</script>")
    assert "<script>alert(1)</script>" not in page


def test_write_report_creates_parent_dirs(tmp_path):
    out = write_report(_stats(), PANEL, tmp_path / "nested" / "r.html")
    assert out.exists() and out.read_text(encoding="utf-8").startswith("<!doctype html>")


def test_cli_demux_report_flag(tmp_path):
    r1, r2 = simulate_reads(PANEL, reads_per_bin=40, target_length=80,
                            error_rate=0.0, paired_end=True, seed=1)
    p1, p2 = tmp_path / "s_R1.fq.gz", tmp_path / "s_R2.fq.gz"
    write_fastq(r1, p1)
    write_fastq(r2, p2)
    cfg = tmp_path / "panel.yaml"
    PANEL.to_yaml(cfg)

    out = tmp_path / "bins"
    rc = main(["demux", "--r1", str(p1), "--r2", str(p2), "--config", str(cfg),
               "--out", str(out), "--report"])
    assert rc == 0
    assert (out / "report.html").exists()

    custom = tmp_path / "custom.html"
    rc = main(["demux", "--r1", str(p1), "--r2", str(p2), "--config", str(cfg),
               "--out", str(tmp_path / "bins2"), "--report", str(custom)])
    assert rc == 0 and custom.exists()

    # without the flag, no report is written
    out3 = tmp_path / "bins3"
    main(["demux", "--r1", str(p1), "--r2", str(p2), "--config", str(cfg), "--out", str(out3)])
    assert not (out3 / "report.html").exists()
