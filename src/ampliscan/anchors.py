"""Locate constant anchor sequences inside reads."""
from __future__ import annotations

from typing import Optional, Tuple

from ampliscan.matcher import hamming


def find_anchor(
    read: str,
    anchor: str,
    max_mismatch: int = 0,
    start: int = 0,
    end: Optional[int] = None,
    from_end: bool = False,
) -> Optional[int]:
    """Sliding-window search for ``anchor`` in ``read[start:end]``.

    Returns the 0-based index of the first occurrence within the mismatch
    budget, or None. If ``from_end`` is True, scan right-to-left and return
    the rightmost occurrence. Mismatch is Hamming (no indels).
    """
    if end is None:
        end = len(read)
    L = len(anchor)
    if L == 0 or end - start < L:
        return None
    read_upper = read.upper()
    anchor_upper = anchor.upper()

    if max_mismatch == 0:
        idx = read_upper.rfind(anchor_upper, start, end) if from_end \
              else read_upper.find(anchor_upper, start, end)
        return idx if idx >= 0 else None

    positions = range(start, end - L + 1)
    if from_end:
        positions = reversed(positions)
    best_pos: Optional[int] = None
    best_dist = max_mismatch + 1
    for i in positions:
        d = hamming(read_upper[i:i + L], anchor_upper)
        if d < best_dist:
            best_dist = d
            best_pos = i
            if d == 0:
                break
    return best_pos if best_dist <= max_mismatch else None


def locate_barcoded_region(
    read: str,
    anchor_5p: str,
    anchor_3p: str,
    barcode_length: int,
    anchor_mismatch: int = 0,
    adapter: str = "",
) -> Optional[Tuple[int, int, str, str, str]]:
    """Find the structure ``ANCHOR5 [F-bc] TARGET [R-bc] ANCHOR3`` in ``read``.

    Returns a tuple ``(fbc_start, rbc_end, fbc, target, rbc)`` where
    ``fbc`` and ``rbc`` are the extracted barcode strings and ``target`` is
    everything between them. Returns None if either anchor is missing or the
    geometry is inconsistent.

    ``adapter``: if given and found after the 5' barcode, the 3' anchor is only
    searched for before it. Reads that run through the insert into the sequencing
    adapter and sample index would otherwise match a stray copy of the anchor
    inside the index (the search takes the rightmost match, which is still the
    rule inside the bounded window, since the insert may contain the anchor).
    """
    pos5 = find_anchor(read, anchor_5p, max_mismatch=anchor_mismatch)
    if pos5 is None:
        return None
    fbc_start = pos5 + len(anchor_5p)
    fbc_end = fbc_start + barcode_length
    if fbc_end > len(read):
        return None

    end3 = None
    if adapter:
        a = read.upper().find(adapter.upper(), fbc_end)
        if a >= 0:
            end3 = a

    pos3 = find_anchor(read, anchor_3p, max_mismatch=anchor_mismatch,
                       start=fbc_end, end=end3, from_end=True)
    if pos3 is None:
        return None
    rbc_end = pos3
    rbc_start = pos3 - barcode_length
    if rbc_start < fbc_end:
        return None

    fbc = read[fbc_start:fbc_end]
    rbc = read[rbc_start:rbc_end]
    target = read[fbc_end:rbc_start]
    return fbc_start, rbc_end + len(anchor_3p), fbc, target, rbc
