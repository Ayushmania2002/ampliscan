"""Parse a draft Excel sheet of amplicons into an ampliscan Panel.

The expected sheet format is what the user's PI sent: a column of "Read_ID"
labels like ``P01_F01.R01`` and a column of "Amplicon" strings with the
literal placeholder ``(X100)`` (or ``(X'100)`` on the reverse strand) where
the variable target sits. The structure scanned for is::

    <adapter> ANCHOR5 <F-barcode> (X100) <R-barcode> ANCHOR3 <adapter>

Anchors are inferred as the 4 nt that immediately precede the F-barcode and
immediately follow the R-barcode in every row. Barcode names are taken from
the Read_ID labels (``F01``, ``R02`` etc.).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from ampliscan.config import Panel, reverse_complement


_TARGET_TOKEN = re.compile(r"\(X'?100\)")


def _parse_amplicon(amplicon: str, barcode_length: int) -> Tuple[str, str, str, str]:
    """Return (anchor5, fbc, rbc, anchor3) for one forward-strand amplicon."""
    m = _TARGET_TOKEN.search(amplicon)
    if m is None:
        raise ValueError(f"Amplicon missing (X100) token: {amplicon!r}")
    before = amplicon[:m.start()]
    after = amplicon[m.end():]
    if len(before) < barcode_length + 1 or len(after) < barcode_length + 1:
        raise ValueError(f"Amplicon too short around (X100): {amplicon!r}")
    fbc = before[-barcode_length:]
    rbc = after[:barcode_length]
    # take 4 nt immediately flanking the barcode as the anchor (configurable)
    anchor_len = 4
    anchor5 = before[-(barcode_length + anchor_len):-barcode_length]
    anchor3 = after[barcode_length:barcode_length + anchor_len]
    return anchor5.upper(), fbc.upper(), rbc.upper(), anchor3.upper()


_LABEL_RE = re.compile(r"(F\d+|R\d+)", re.IGNORECASE)


def _extract_label_pair(read_id: str) -> Tuple[Optional[str], Optional[str]]:
    """Pull F-name and R-name out of an ID like 'P01_F01.R01' or 'P01_R01.F01'."""
    tokens = _LABEL_RE.findall(read_id)
    f_name = next((t.upper() for t in tokens if t.upper().startswith("F")), None)
    r_name = next((t.upper() for t in tokens if t.upper().startswith("R")), None)
    return f_name, r_name


def panel_from_excel(
    path: Union[str, Path],
    sheet_name: Union[str, int] = 0,
    barcode_length: int = 4,
    anchor_length: int = 4,
    name: str = "panel_from_excel",
) -> Panel:
    """Build a Panel by scanning an Excel amplicon table.

    Only forward-strand rows (those containing ``(X100)``) are used to define
    the panel; reverse-strand rows are consistency-checked if present.
    """
    import pandas as pd

    df = pd.read_excel(path, sheet_name=sheet_name, header=0)
    if df.shape[1] < 2:
        raise ValueError("Excel sheet must have at least 2 columns: Read_ID, Amplicon")
    id_col, amp_col = df.columns[0], df.columns[1]

    forward_barcodes: Dict[str, str] = {}
    reverse_barcodes: Dict[str, str] = {}
    anchor5_set: set[str] = set()
    anchor3_set: set[str] = set()

    for _, row in df.iterrows():
        rid = str(row[id_col]).strip().rstrip("'")
        amp = str(row[amp_col]).strip()
        if "(X100)" not in amp:
            continue  # skip reverse-strand rows
        a5, fbc, rbc, a3 = _parse_amplicon(amp, barcode_length)
        anchor5_set.add(a5)
        anchor3_set.add(a3)
        f_name, r_name = _extract_label_pair(rid)
        if f_name:
            forward_barcodes.setdefault(f_name, fbc)
            if forward_barcodes[f_name] != fbc:
                raise ValueError(
                    f"Inconsistent F barcode for {f_name}: "
                    f"{forward_barcodes[f_name]} vs {fbc}"
                )
        if r_name:
            reverse_barcodes.setdefault(r_name, rbc)
            if reverse_barcodes[r_name] != rbc:
                raise ValueError(
                    f"Inconsistent R barcode for {r_name}: "
                    f"{reverse_barcodes[r_name]} vs {rbc}"
                )

    if len(anchor5_set) != 1:
        raise ValueError(f"Multiple 5' anchors detected: {anchor5_set}")
    if len(anchor3_set) != 1:
        raise ValueError(f"Multiple 3' anchors detected: {anchor3_set}")
    if not forward_barcodes:
        raise ValueError("No F barcodes found in Excel")
    if not reverse_barcodes:
        raise ValueError("No R barcodes found in Excel")

    return Panel(
        name=name,
        forward_5p_anchor=next(iter(anchor5_set)),
        forward_3p_anchor=next(iter(anchor3_set)),
        forward_barcodes=dict(sorted(forward_barcodes.items())),
        reverse_barcodes=dict(sorted(reverse_barcodes.items())),
        barcode_length=barcode_length,
        notes=f"parsed from {Path(path).name}",
    )
