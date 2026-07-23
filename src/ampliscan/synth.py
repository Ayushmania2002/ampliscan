"""Synthetic read generator for testing/demo."""
from __future__ import annotations

import random
from typing import Iterator, List, Optional, Tuple

from ampliscan.config import Panel, reverse_complement
from ampliscan.io import Record


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=n))


def mutate(seq: str, error_rate: float, rng: random.Random) -> str:
    """Per-base substitution noise (simple model — no indels)."""
    if error_rate <= 0:
        return seq
    out = []
    bases = "ACGT"
    for b in seq:
        if rng.random() < error_rate:
            alt = bases.replace(b.upper(), "")
            out.append(rng.choice(alt))
        else:
            out.append(b)
    return "".join(out)


def simulate_reads(
    panel: Panel,
    reads_per_bin: int = 100,
    target_length: int = 100,
    error_rate: float = 0.005,
    paired_end: bool = True,
    include_adapter_flanks: bool = True,
    seed: Optional[int] = 42,
) -> Tuple[List[Record], Optional[List[Record]]]:
    """Generate synthetic reads for every F x R bin in the panel.

    Returns (r1_records, r2_records). r2 is None if not paired_end.
    """
    rng = random.Random(seed)
    p5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" if include_adapter_flanks else ""
    p7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGCGATCTCGTATGCCGTCTTCTGCTTG" if include_adapter_flanks else ""

    r1_list: List[Record] = []
    r2_list: List[Record] = []
    idx = 0
    for f_name, f_seq in panel.forward_barcodes.items():
        for r_name, r_seq in panel.reverse_barcodes.items():
            for k in range(reads_per_bin):
                target = random_dna(target_length, rng)
                molecule = (
                    p5
                    + panel.forward_5p_anchor
                    + f_seq
                    + target
                    + r_seq
                    + panel.forward_3p_anchor
                    + p7
                )
                r1_seq = mutate(molecule, error_rate, rng)
                rid = f"sim{idx:06d}_{f_name}_{r_name}_{k}"
                r1_list.append(Record(id=rid, seq=r1_seq, qual="I" * len(r1_seq)))
                if paired_end:
                    r2_seq = mutate(reverse_complement(molecule), error_rate, rng)
                    r2_list.append(Record(id=rid, seq=r2_seq, qual="I" * len(r2_seq)))
                idx += 1
    rng.shuffle(r1_list)
    if paired_end:
        # keep r2 paired with r1
        id_to_r2 = {r.id: r for r in r2_list}
        r2_list = [id_to_r2[r.id] for r in r1_list]
        return r1_list, r2_list
    return r1_list, None
