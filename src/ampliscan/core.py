"""Demultiplexing pipeline: read in, bin out."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Iterator, Optional, Tuple

from ampliscan.config import Panel, MatchPolicy, reverse_complement
from ampliscan.anchors import locate_barcoded_region
from ampliscan.matcher import BarcodeMatcher
from ampliscan.io import Record


REASON_NO_ANCHOR = "no_anchor"
REASON_BC_UNKNOWN = "barcode_unmatched"
REASON_BC_AMBIGUOUS = "barcode_ambiguous"
REASON_PAIR_MISMATCH = "paired_end_disagreement"
REASON_LOW_QUALITY = "low_quality"


@dataclass
class DemuxResult:
    """Outcome of demultiplexing one read or read pair."""
    read_id: str
    bin_name: Optional[str]          # "F01_R01" if assigned, else None
    forward_bc: Optional[str]        # name, e.g. "F01"
    reverse_bc: Optional[str]        # name, e.g. "R01"
    target_seq: Optional[str]        # extracted insert
    target_qual: Optional[str] = None
    strand: Optional[str] = None     # "+" or "-" (orientation of the returned/winning result)
    reason: Optional[str] = None     # failure reason if unassigned
    fbc_distance: int = 0
    rbc_distance: int = 0
    # Independent per-mate bin assignment for paired-end input (None for
    # single-end). Unlike `strand`, which only reflects whichever mate's
    # result was returned, these are populated from *each* mate's own
    # demultiplexing regardless of pairing outcome -- e.g. for a QC heatmap
    # comparing what R1 alone vs R2 alone would assign.
    r1_bin_name: Optional[str] = None
    r2_bin_name: Optional[str] = None


def _try_orientation(
    seq: str,
    qual: Optional[str],
    panel: Panel,
    fwd_matcher: BarcodeMatcher,
    rev_matcher: BarcodeMatcher,
    policy: MatchPolicy,
    strand: str,
) -> Optional[DemuxResult]:
    """Try to demux ``seq`` as if it were on ``strand``. Returns None on no anchor."""
    if strand == "+":
        a5, a3 = panel.forward_5p_anchor, panel.forward_3p_anchor
    else:
        a5, a3 = panel.reverse_5p_anchor, panel.reverse_3p_anchor

    hit = locate_barcoded_region(
        seq, a5, a3, panel.barcode_length, anchor_mismatch=policy.anchor_mismatch
    )
    if hit is None:
        return None

    _, _, fbc_seen, target, rbc_seen = hit
    if strand == "+":
        fbc_query, rbc_query = fbc_seen, rbc_seen
    else:
        # On the reverse-strand read, the barcode adjacent to ANCHOR5(rev) is
        # actually the RC of the R-side barcode, and the barcode adjacent to
        # ANCHOR3(rev) is the RC of the F-side barcode.
        fbc_query = reverse_complement(rbc_seen)
        rbc_query = reverse_complement(fbc_seen)
        target = reverse_complement(target)

    f_name, f_d, f_amb = fwd_matcher.match(fbc_query)
    r_name, r_d, r_amb = rev_matcher.match(rbc_query)

    if f_amb or r_amb:
        reason = REASON_BC_AMBIGUOUS
        bin_name = None
    elif f_name is None or r_name is None:
        reason = REASON_BC_UNKNOWN
        bin_name = None
    else:
        reason = None
        bin_name = f"{f_name}_{r_name}"

    return DemuxResult(
        read_id="",  # filled in by caller
        bin_name=bin_name,
        forward_bc=f_name,
        reverse_bc=r_name,
        target_seq=target,
        target_qual=None,  # quality slicing for now is best-effort, skipped in v0.1
        strand=strand,
        reason=reason,
        fbc_distance=f_d,
        rbc_distance=r_d,
    )


def demux_record(
    record: Record,
    panel: Panel,
    policy: Optional[MatchPolicy] = None,
    mate: Optional[Record] = None,
) -> DemuxResult:
    """Demultiplex one record (single-end) or one read pair (R1 + mate=R2)."""
    policy = policy or MatchPolicy()
    fwd_matcher = BarcodeMatcher(
        panel.forward_barcodes, policy.barcode_mismatch, policy.allow_indels
    )
    rev_matcher = BarcodeMatcher(
        panel.reverse_barcodes, policy.barcode_mismatch, policy.allow_indels
    )
    return _demux_with_matchers(record, panel, policy, fwd_matcher, rev_matcher, mate)


def _demux_with_matchers(record, panel, policy, fwd_matcher, rev_matcher, mate):
    r1 = _demux_single(record, panel, policy, fwd_matcher, rev_matcher)
    if mate is None:
        r1.r1_bin_name = r1.bin_name
        return r1

    r2 = _demux_single(mate, panel, policy, fwd_matcher, rev_matcher)

    if r1.bin_name and r2.bin_name:
        if r1.bin_name == r2.bin_name:
            r1.read_id = record.id
            r1.r1_bin_name = r1.bin_name
            r1.r2_bin_name = r2.bin_name
            return r1
        # disagreement -- mark unassigned
        return DemuxResult(
            read_id=record.id,
            bin_name=None,
            forward_bc=None,
            reverse_bc=None,
            target_seq=r1.target_seq,
            strand=r1.strand,
            reason=REASON_PAIR_MISMATCH,
            fbc_distance=r1.fbc_distance,
            rbc_distance=r1.rbc_distance,
            r1_bin_name=r1.bin_name,
            r2_bin_name=r2.bin_name,
        )
    # one side worked, use it
    winner = r1 if r1.bin_name else r2
    winner.read_id = record.id
    winner.r1_bin_name = r1.bin_name
    winner.r2_bin_name = r2.bin_name
    return winner


def _demux_single(record, panel, policy, fwd_matcher, rev_matcher):
    for strand in ("+", "-"):
        res = _try_orientation(
            record.seq, record.qual, panel, fwd_matcher, rev_matcher, policy, strand
        )
        if res is not None:
            res.read_id = record.id
            if res.bin_name is not None:
                return res
            best = res  # remember failure to report
    if 'best' in locals():
        return best
    return DemuxResult(
        read_id=record.id,
        bin_name=None,
        forward_bc=None,
        reverse_bc=None,
        target_seq=None,
        reason=REASON_NO_ANCHOR,
    )


def demux_records(
    records: Iterable[Record],
    panel: Panel,
    policy: Optional[MatchPolicy] = None,
    mates: Optional[Iterable[Record]] = None,
) -> Iterator[DemuxResult]:
    """Stream a demux over an iterable of records (optionally paired)."""
    policy = policy or MatchPolicy()
    fwd_matcher = BarcodeMatcher(
        panel.forward_barcodes, policy.barcode_mismatch, policy.allow_indels
    )
    rev_matcher = BarcodeMatcher(
        panel.reverse_barcodes, policy.barcode_mismatch, policy.allow_indels
    )
    if mates is None:
        for rec in records:
            yield _demux_with_matchers(rec, panel, policy, fwd_matcher, rev_matcher, None)
    else:
        for r1, r2 in zip(records, mates):
            yield _demux_with_matchers(r1, panel, policy, fwd_matcher, rev_matcher, r2)
