"""Approximate string matching for barcodes and anchors."""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple


def hamming(a: str, b: str) -> int:
    """Hamming distance; returns len(a) if lengths differ."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def levenshtein(a: str, b: str, max_dist: Optional[int] = None) -> int:
    """Banded Levenshtein with early exit at max_dist."""
    if a == b:
        return 0
    if max_dist is not None and abs(len(a) - len(b)) > max_dist:
        return max_dist + 1
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        curr = [i] + [0] * len(b)
        row_min = curr[0]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
            if curr[j] < row_min:
                row_min = curr[j]
        if max_dist is not None and row_min > max_dist:
            return max_dist + 1
        prev = curr
    return prev[-1]


class BarcodeMatcher:
    """Look up a query barcode against a panel with bounded mismatches.

    The matcher is unambiguous only when exactly one panel barcode lies within
    the mismatch budget. Otherwise the result is ``None`` (no hit) or a
    tie/ambiguous result (multiple equally-close hits).
    """

    def __init__(
        self,
        barcodes: Dict[str, str],
        max_mismatch: int = 1,
        allow_indels: bool = False,
    ):
        self.barcodes = dict(barcodes)
        self.max_mismatch = max_mismatch
        self.allow_indels = allow_indels
        # exact-lookup fast path
        self._exact = {seq: name for name, seq in self.barcodes.items()}

    def match(self, query: str) -> Tuple[Optional[str], int, bool]:
        """Return (best_name, distance, ambiguous).

        ``best_name`` is None if no barcode is within the mismatch budget.
        ``ambiguous`` is True if two or more barcodes tie at the best distance.
        """
        query = query.upper()
        if query in self._exact:
            return self._exact[query], 0, False

        distfn = levenshtein if self.allow_indels else hamming
        best_name: Optional[str] = None
        best_dist = self.max_mismatch + 1
        ambiguous = False
        for name, seq in self.barcodes.items():
            d = distfn(query, seq) if not self.allow_indels else distfn(query, seq, self.max_mismatch)
            if d < best_dist:
                best_dist = d
                best_name = name
                ambiguous = False
            elif d == best_dist and best_name is not None:
                ambiguous = True

        if best_name is None or best_dist > self.max_mismatch:
            return None, best_dist, False
        return best_name, best_dist, ambiguous
