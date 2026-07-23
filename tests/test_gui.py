"""Tests for GUI-support helpers (no display required).

The full Tk window can't be exercised on a headless CI runner, so these cover
the pure helper functions the GUI relies on. A display-dependent smoke test is
skipped automatically when Tk can't open a window.
"""
import pytest

from ampliscan.gui import _parse_barcodes, _viridis, _require_tk


def test_parse_barcodes_separators():
    assert _parse_barcodes("F01\tGCGT\nF02\tGTAG") == {"F01": "GCGT", "F02": "GTAG"}
    assert _parse_barcodes("F01,GCGT\nF02,GTAG") == {"F01": "GCGT", "F02": "GTAG"}
    assert _parse_barcodes("F01 GCGT") == {"F01": "GCGT"}


def test_parse_barcodes_uppercases_and_skips_blanks():
    assert _parse_barcodes("F01 gcgt\n\n  \nF02 gtag") == {"F01": "GCGT", "F02": "GTAG"}


def test_viridis_endpoints():
    assert _viridis(0.0) == "#440154"
    assert _viridis(1.0) == "#fde725"
    # clamps out-of-range input
    assert _viridis(-5) == "#440154"
    assert _viridis(5) == "#fde725"


def test_gui_constructs_if_display_available():
    if not _require_tk():
        pytest.skip("tkinter not available in this build")
    try:
        import tkinter as tk
        from ampliscan.gui import AmpliscanGUI
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available")
    root.withdraw()
    AmpliscanGUI(root)          # must not raise
    root.update()               # force a render pass
    root.destroy()


def test_gui_has_two_heatmap_canvases_and_renders_results():
    if not _require_tk():
        pytest.skip("tkinter not available in this build")
    try:
        import tkinter as tk
        from ampliscan.gui import AmpliscanGUI
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available")
    root.withdraw()
    app = AmpliscanGUI(root)
    root.update()

    # two distinct canvases power the R1/R2 comparison heatmaps
    assert app.canvas_r1 is not app.canvas_r2

    from ampliscan.parallel import DemuxStats
    from collections import Counter
    stats = DemuxStats(total=10, assigned=8)
    stats.bin_counts_r1 = Counter({"F01_R01": 5, "F02_R02": 3})
    stats.bin_counts_r2 = Counter({"F01_R01": 4, "F02_R02": 4})
    app._show_results(stats)   # must not raise
    root.update()

    assert app._last_bins_r1 == stats.bin_counts_r1
    assert app._last_bins_r2 == stats.bin_counts_r2
    root.destroy()
