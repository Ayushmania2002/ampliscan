"""Self-contained HTML report for a demultiplexing run.

Standard library only. The output is a single HTML file (inline CSS + SVG, no
external requests) containing the F x R heatmap, per-mate (R1 vs R2) heatmaps
for paired-end runs, a donut chart of unassigned reasons, and per-bin tables.
"""
from __future__ import annotations

import html
import math
from datetime import date
from pathlib import Path
from typing import Iterable, Optional, Union

from ampliscan import __version__

_VIRIDIS = [(68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98), (253, 231, 37)]
_PALETTE = ["#2563eb", "#d97706", "#0d9488", "#9333ea", "#dc2626", "#65a30d"]
_REASON_LABELS = {
    "no_anchor": "No anchor found",
    "barcode_unmatched": "Barcode unmatched",
    "barcode_ambiguous": "Barcode ambiguous",
    "paired_end_disagreement": "R1/R2 disagree",
}


def _viridis(t: float):
    t = max(0.0, min(1.0, t))
    s = t * (len(_VIRIDIS) - 1)
    i = min(int(s), len(_VIRIDIS) - 2)
    f = s - i
    a, b = _VIRIDIS[i], _VIRIDIS[i + 1]
    return tuple(round(a[k] + (b[k] - a[k]) * f) for k in range(3))


def _fmt(n: int) -> str:
    return f"{n:,}"


def _heatmap(counts, f_names, r_names, vmax: int, denom: int, compact: bool = False) -> str:
    esc = html.escape
    rows = ["<tr><th></th>" + "".join(f'<th scope="col">{esc(r)}</th>' for r in r_names) + "</tr>"]
    for f in f_names:
        cells = [f'<th scope="row">{esc(f)}</th>']
        for r in r_names:
            v = counts.get(f"{f}_{r}", 0)
            rgb = _viridis(v / vmax if vmax else 0)
            lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
            fg = "#111" if lum > 140 else "#fff"
            pct = 100 * v / denom if denom else 0
            small = "" if compact else f"<small>{pct:.1f}%</small>"
            cells.append(
                f'<td style="background:rgb{rgb};color:{fg}" '
                f'title="{esc(f)}_{esc(r)}: {_fmt(v)} reads ({pct:.2f}%)">'
                f"<b>{_fmt(v)}</b>{small}</td>"
            )
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return '<table class="heat">' + "".join(rows) + "</table>"


def _colorbar(vmax: int) -> str:
    stops = ", ".join(f"rgb{_viridis(i / 10)}" for i in range(11))
    return (f'<div class="cbar"><span>0</span>'
            f'<div style="background:linear-gradient(90deg,{stops})"></div>'
            f"<span>{_fmt(vmax)} reads</span></div>")


def _donut(reasons, unassigned: int):
    items = sorted(reasons.items(), key=lambda kv: -kv[1])
    cx = cy = 110
    R, r = 100, 58
    angle = -math.pi / 2
    paths, legend = [], []
    for i, (k, v) in enumerate(items):
        frac = v / unassigned if unassigned else 0
        a2 = angle + 2 * math.pi * frac
        if frac >= 0.9999:  # a single full-circle arc degenerates; nudge it
            a2 = angle + 2 * math.pi * 0.9999
        large = 1 if frac > 0.5 else 0
        x1, y1 = cx + R * math.cos(angle), cy + R * math.sin(angle)
        x2, y2 = cx + R * math.cos(a2), cy + R * math.sin(a2)
        x3, y3 = cx + r * math.cos(a2), cy + r * math.sin(a2)
        x4, y4 = cx + r * math.cos(angle), cy + r * math.sin(angle)
        col = _PALETTE[i % len(_PALETTE)]
        label = _REASON_LABELS.get(k, k)
        paths.append(
            f'<path d="M{x1:.2f},{y1:.2f} A{R},{R} 0 {large} 1 {x2:.2f},{y2:.2f} '
            f'L{x3:.2f},{y3:.2f} A{r},{r} 0 {large} 0 {x4:.2f},{y4:.2f} Z" '
            f'fill="{col}" stroke="#fff" stroke-width="2">'
            f"<title>{html.escape(label)}: {_fmt(v)} ({100 * frac:.1f}% of unassigned)</title></path>"
        )
        legend.append(
            f'<tr><td><span class="sw" style="background:{col}"></span>{html.escape(label)}<br>'
            f'<code>{html.escape(k)}</code></td><td class="num">{_fmt(v)}</td>'
            f'<td class="num">{100 * frac:.1f}%</td></tr>'
        )
        angle = a2
    svg = (f'<svg viewBox="0 0 220 220" role="img" aria-label="Unassigned reads by reason">'
           f'{"".join(paths)}<text x="110" y="106" text-anchor="middle" class="dn">{_fmt(unassigned)}</text>'
           f'<text x="110" y="126" text-anchor="middle" class="dl">unassigned</text></svg>')
    return svg, "".join(legend)


_CSS = """
:root{--bg:#f6f7f9;--card:#fff;--fg:#111827;--mut:#6b7280;--bd:#e5e7eb}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 "Segoe UI",system-ui,sans-serif}
header{background:var(--card);border-bottom:1px solid var(--bd);padding:22px 20px}
.wrap{max-width:1040px;margin:0 auto;padding:0 20px}
header .wrap{padding:0}
h1{margin:0;font-size:24px}
h2{font-size:13px;letter-spacing:.06em;text-transform:uppercase;color:var(--mut);margin:0 0 14px}
.sub{color:var(--mut);margin-top:2px;font-size:14px}
main{padding:22px 0 40px}
.card{background:var(--card);border:1px solid var(--bd);border-radius:10px;padding:20px;margin-bottom:18px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:14px}
.tile{background:var(--bg);border-radius:8px;padding:14px 16px}
.tile b{display:block;font-size:26px}.tile span{color:var(--mut);font-size:13px}
.good b{color:#059669}
table{border-collapse:collapse;width:100%}
.heat{table-layout:fixed}
.heat th,.heat td{text-align:center}
.heat th{color:var(--mut);font-size:12px;padding:6px 4px;font-weight:600}
.heat th:first-child{width:48px}
.heat td{padding:12px 4px;border:2px solid #fff;border-radius:4px}
.heat td b{display:block;font-size:15px}.heat td small{opacity:.85;font-size:11px}
.pair .heat td{padding:9px 2px}.pair .heat td b{font-size:12.5px}
.cbar{display:flex;align-items:center;gap:10px;margin-top:12px;font-size:12px;color:var(--mut)}
.cbar div{flex:1;height:10px;border-radius:5px;max-width:320px}
.two{display:grid;grid-template-columns:240px 1fr;gap:26px;align-items:center}
.two svg{width:100%;height:auto}
.dn{font-size:21px;font-weight:700;fill:var(--fg)}.dl{font-size:11px;fill:var(--mut)}
.sw{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:8px}
td,th{padding:7px 8px;text-align:left;border-bottom:1px solid var(--bd);font-size:14px}
.num{text-align:right;font-variant-numeric:tabular-nums}
code{background:var(--bg);padding:1px 5px;border-radius:4px;font-size:12px;color:var(--mut)}
.cols,.pair{display:grid;grid-template-columns:1fr 1fr;gap:18px}
.pair h3{margin:0 0 8px;font-size:13px;color:var(--mut);font-weight:600}
.note{font-size:13.5px;color:var(--mut)}.note li{margin:4px 0}
footer{text-align:center;color:var(--mut);font-size:12.5px;padding-bottom:30px}
@media(max-width:760px){.two,.cols,.pair{grid-template-columns:1fr}}
@media print{body{background:#fff}.card{break-inside:avoid}}
"""


def build_report_html(stats, panel, title: str = "", params: Optional[dict] = None) -> str:
    """Return the report as an HTML string."""
    esc = html.escape
    bins = dict(stats.bin_counts)
    total, assigned, unassigned = stats.total, stats.assigned, stats.unassigned
    f_names = list(panel.forward_barcodes)
    r_names = list(panel.reverse_barcodes)
    vmax = max(bins.values(), default=0)

    main_heat = _heatmap(bins, f_names, r_names, vmax, assigned) + _colorbar(vmax)

    # Per-mate heatmaps only make sense for paired-end runs.
    c1, c2 = dict(stats.bin_counts_r1), dict(stats.bin_counts_r2)
    mate_card = ""
    if c1 and c2:
        mmax = max(max(c1.values()), max(c2.values()))
        s1, s2 = sum(c1.values()), sum(c2.values())
        mate_card = (
            '<div class="card"><h2>R1 vs R2: independent assignment</h2><div class="pair">'
            f'<div><h3>R1 alone ({_fmt(s1)} reads assigned)</h3>{_heatmap(c1, f_names, r_names, mmax, s1, True)}</div>'
            f'<div><h3>R2 alone ({_fmt(s2)} reads assigned)</h3>{_heatmap(c2, f_names, r_names, mmax, s2, True)}</div>'
            f"</div>{_colorbar(mmax)}"
            '<p class="note">Each mate is demultiplexed on its own, before the pair is reconciled. '
            "Both panels use the same colour scale, so similar patterns mean the mates agree; "
            "a bin that is much darker in one panel points to a problem with that mate's anchor or barcode.</p></div>"
        )

    donut, legend = _donut(stats.reason_counts, unassigned) if unassigned else ("", "")
    if unassigned:
        reasons_card = (
            '<div class="card"><h2>Why reads were not assigned</h2>'
            f'<div class="two">{donut}<table><thead><tr><th>Reason</th><th class="num">Reads</th>'
            f'<th class="num">Share</th></tr></thead><tbody>{legend}</tbody></table></div>'
            '<ul class="note">'
            f"<li><b>No anchor found</b>: neither flanking anchor (<code>{esc(panel.forward_5p_anchor)}</code> / "
            f"<code>{esc(panel.forward_3p_anchor)}</code>) was located.</li>"
            "<li><b>Barcode unmatched</b>: anchors found, but the barcode is further than the mismatch budget "
            "from every panel barcode.</li>"
            "<li><b>Barcode ambiguous</b>: the barcode is equally close to two panel barcodes, so no guess was made.</li>"
            "<li><b>R1/R2 disagree</b>: the two mates were assigned to different samples; the pair is withheld.</li>"
            "</ul></div>"
        )
    else:
        reasons_card = '<div class="card"><h2>Why reads were not assigned</h2><p class="note">Every read was assigned.</p></div>'

    bin_rows = "".join(
        f'<tr><td>{esc(k)}</td><td class="num">{_fmt(v)}</td>'
        f'<td class="num">{100 * v / assigned if assigned else 0:.2f}%</td>'
        f'<td class="num">{100 * v / total if total else 0:.2f}%</td></tr>'
        for k, v in sorted(bins.items(), key=lambda kv: -kv[1])
    )
    vals = sorted(bins.values())
    fold_note = (f'<p class="note">Highest bin is {vals[-1] / max(vals[0], 1):.1f}&times; the lowest '
                 f"({_fmt(vals[-1])} vs {_fmt(vals[0])}).</p>") if vals else ""

    panel_rows = "".join(
        f"<tr><td>{esc(k)}</td><td><code>{esc(v)}</code></td></tr>"
        for k, v in {**panel.forward_barcodes, **panel.reverse_barcodes}.items()
    )
    param_note = ""
    if params:
        param_note = "<p class=\"note\">" + " &middot; ".join(
            f"{esc(str(k))}: <code>{esc(str(v))}</code>" for k, v in params.items()) + "</p>"

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ampliscan report: {esc(title or panel.name)}</title>
<style>{_CSS}</style></head><body>
<header><div class="wrap">
<h1>ampliscan demultiplexing report</h1>
<div class="sub">{esc(title)}{' &middot; ' if title else ''}panel <code>{esc(panel.name)}</code> &middot; generated {date.today().isoformat()}</div>
</div></header>
<main class="wrap">
<div class="card"><h2>Summary</h2><div class="tiles">
<div class="tile"><b>{_fmt(total)}</b><span>reads / read pairs processed</span></div>
<div class="tile good"><b>{100 * stats.assigned_fraction:.1f}%</b><span>assigned to a sample ({_fmt(assigned)})</span></div>
<div class="tile"><b>{_fmt(unassigned)}</b><span>unassigned ({100 * unassigned / total if total else 0:.1f}%)</span></div>
<div class="tile"><b>{len(bins)}</b><span>of {len(f_names) * len(r_names)} sample bins populated</span></div>
</div></div>
<div class="card"><h2>Reads per forward &times; reverse barcode bin</h2>{main_heat}
<p class="note">Cell colour is scaled to the largest bin. Percentages are of assigned reads. Hover a cell for its exact count.</p></div>
{mate_card}
{reasons_card}
<div class="cols">
<div class="card"><h2>Reads per bin</h2><table><thead><tr><th>Bin</th><th class="num">Reads</th>
<th class="num">% assigned</th><th class="num">% total</th></tr></thead><tbody>{bin_rows}</tbody></table>{fold_note}</div>
<div class="card"><h2>Panel used</h2><table><thead><tr><th>Name</th><th>Barcode</th></tr></thead><tbody>{panel_rows}</tbody></table>
<p class="note">Anchors: 5&prime; <code>{esc(panel.forward_5p_anchor)}</code>, 3&prime; <code>{esc(panel.forward_3p_anchor)}</code>
&middot; barcode length {panel.barcode_length}.</p>{param_note}</div>
</div></main>
<footer>Generated with ampliscan {esc(__version__)} &middot; doi:10.5281/zenodo.21505581</footer>
</body></html>
"""


def write_report(stats, panel, path: Union[str, Path], title: str = "",
                 params: Optional[dict] = None) -> Path:
    """Write the HTML report to ``path`` and return the resolved path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_report_html(stats, panel, title=title, params=params), encoding="utf-8")
    return path.resolve()
