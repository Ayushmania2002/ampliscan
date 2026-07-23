"""Panel configuration: user-defined anchors and barcodes."""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, Optional, Union, List, Tuple


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class MatchPolicy:
    """How tolerant the matcher should be.

    anchor_mismatch: substitutions allowed in the constant anchor (default 0).
    barcode_mismatch: substitutions allowed in the variable barcode (default 1).
    allow_indels: if True, use Levenshtein distance; otherwise Hamming.
    min_anchor_quality: minimum mean Phred over the anchor region (0 = no filter).
    min_barcode_quality: minimum mean Phred over the barcode region (0 = no filter).
    """
    anchor_mismatch: int = 0
    barcode_mismatch: int = 1
    allow_indels: bool = False
    min_anchor_quality: int = 0
    min_barcode_quality: int = 0


@dataclass
class Panel:
    """A demultiplexing panel.

    Everything is user-defined: anchor sequences, barcode length, and the
    barcode dictionaries themselves. Reverse-side anchors default to the
    reverse complement of the forward-side anchors but can be overridden if
    the library design is asymmetric.
    """
    name: str
    forward_5p_anchor: str
    forward_3p_anchor: str
    forward_barcodes: Dict[str, str]
    reverse_barcodes: Dict[str, str]
    barcode_length: Optional[int] = None  # inferred if None
    reverse_5p_anchor: Optional[str] = None  # = revcomp(forward_3p) if None
    reverse_3p_anchor: Optional[str] = None  # = revcomp(forward_5p) if None
    target_length_range: Optional[Tuple[int, int]] = None
    notes: str = ""

    def __post_init__(self) -> None:
        self.forward_5p_anchor = self.forward_5p_anchor.upper()
        self.forward_3p_anchor = self.forward_3p_anchor.upper()
        self.forward_barcodes = {k: v.upper() for k, v in self.forward_barcodes.items()}
        self.reverse_barcodes = {k: v.upper() for k, v in self.reverse_barcodes.items()}

        if self.reverse_5p_anchor is None:
            self.reverse_5p_anchor = reverse_complement(self.forward_3p_anchor)
        else:
            self.reverse_5p_anchor = self.reverse_5p_anchor.upper()
        if self.reverse_3p_anchor is None:
            self.reverse_3p_anchor = reverse_complement(self.forward_5p_anchor)
        else:
            self.reverse_3p_anchor = self.reverse_3p_anchor.upper()

        lengths = {len(v) for v in self.forward_barcodes.values()} | \
                  {len(v) for v in self.reverse_barcodes.values()}
        if not lengths:
            raise ValueError("Panel has no barcodes")
        if len(lengths) > 1 and self.barcode_length is None:
            raise ValueError(
                f"Barcodes have inconsistent lengths {sorted(lengths)}; "
                "set barcode_length explicitly or use uniform lengths"
            )
        if self.barcode_length is None:
            self.barcode_length = lengths.pop()

        for name, seq in {**self.forward_barcodes, **self.reverse_barcodes}.items():
            if len(seq) != self.barcode_length:
                raise ValueError(
                    f"Barcode {name!r} has length {len(seq)} but panel "
                    f"barcode_length is {self.barcode_length}"
                )
            if set(seq) - set("ACGTN"):
                raise ValueError(f"Barcode {name!r} has non-DNA characters: {seq!r}")

        if len(set(self.forward_barcodes.values())) != len(self.forward_barcodes):
            raise ValueError("Duplicate forward barcode sequences")
        if len(set(self.reverse_barcodes.values())) != len(self.reverse_barcodes):
            raise ValueError("Duplicate reverse barcode sequences")

    @property
    def bin_names(self) -> List[str]:
        """All F x R combinations as bin identifiers."""
        return [f"{f}_{r}" for f in self.forward_barcodes for r in self.reverse_barcodes]

    @classmethod
    def from_dict(cls, d: dict) -> "Panel":
        tlr = d.get("target_length_range")
        if tlr is not None:
            tlr = tuple(tlr)
        return cls(
            name=d.get("name", "unnamed"),
            forward_5p_anchor=d["forward_5p_anchor"],
            forward_3p_anchor=d["forward_3p_anchor"],
            forward_barcodes=dict(d["forward_barcodes"]),
            reverse_barcodes=dict(d["reverse_barcodes"]),
            barcode_length=d.get("barcode_length"),
            reverse_5p_anchor=d.get("reverse_5p_anchor"),
            reverse_3p_anchor=d.get("reverse_3p_anchor"),
            target_length_range=tlr,
            notes=d.get("notes", ""),
        )

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "Panel":
        import yaml
        with open(path) as fh:
            return cls.from_dict(yaml.safe_load(fh))

    def to_yaml(self, path: Union[str, Path]) -> None:
        import yaml
        data = asdict(self)
        if data["target_length_range"] is not None:
            data["target_length_range"] = list(data["target_length_range"])
        with open(path, "w") as fh:
            yaml.safe_dump(data, fh, sort_keys=False)

    def describe(self) -> str:
        lines = [
            f"Panel: {self.name}",
            f"  5' anchor (fwd): {self.forward_5p_anchor}",
            f"  3' anchor (fwd): {self.forward_3p_anchor}",
            f"  5' anchor (rev): {self.reverse_5p_anchor}",
            f"  3' anchor (rev): {self.reverse_3p_anchor}",
            f"  barcode length:  {self.barcode_length}",
            f"  forward barcodes ({len(self.forward_barcodes)}): "
            + ", ".join(f"{k}={v}" for k, v in self.forward_barcodes.items()),
            f"  reverse barcodes ({len(self.reverse_barcodes)}): "
            + ", ".join(f"{k}={v}" for k, v in self.reverse_barcodes.items()),
            f"  total bins: {len(self.forward_barcodes) * len(self.reverse_barcodes)}",
        ]
        return "\n".join(lines)
