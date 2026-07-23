"""ampliscan: anchor-based amplicon demultiplexer."""

from ampliscan.config import Panel, MatchPolicy
from ampliscan.core import demux_record, demux_records, DemuxResult
from ampliscan.io import iter_records, Record, write_bins
from ampliscan.synth import simulate_reads
from ampliscan.excel_panel import panel_from_excel
from ampliscan.parallel import demux_to_bins, DemuxStats

__version__ = "1.1.0"
__all__ = [
    "Panel",
    "MatchPolicy",
    "Record",
    "DemuxResult",
    "DemuxStats",
    "demux_record",
    "demux_records",
    "demux_to_bins",
    "iter_records",
    "write_bins",
    "simulate_reads",
    "panel_from_excel",
]
