<p align="center">
  <img src="assets/cover.png" alt="ampliscan — sort sequencing reads into per-sample bins by their barcodes" width="820">
</p>

<h1 align="center">ampliscan</h1>

<p align="center">
  <b>Anchor-based amplicon demultiplexer — sort sequencing reads into per-sample bins by their barcodes.</b>
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/python-%E2%89%A53.9-blue">
  <img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green">
  <img alt="Interfaces" src="https://img.shields.io/badge/interfaces-CLI%20%2B%20GUI-8a2be2">
  <!-- After you publish on Zenodo, replace XXXXXXX and uncomment: -->
  <!-- <a href="https://doi.org/10.5281/zenodo.XXXXXXX"><img alt="DOI" src="https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg"></a> -->
</p>

---

ampliscan splits pooled sequencing reads into per-sample files by matching short
barcodes that sit immediately next to user-defined constant **anchor** sequences
(e.g. `gctt` and `acag`). Using the anchor + barcode together — an
*octanucleotide* in the common 4+4 case — gives much better discrimination than
the barcode alone, because the constant anchor pins the alignment even when the
barcode region has a sequencing error.

```
[adapter]  ANCHOR5  [F-barcode]  TARGET  [R-barcode]  ANCHOR3  [adapter]
           └──── octanucleotide ────┘              └──── octanucleotide ────┘
```

It's built for the second, **inner** level of demultiplexing that your
sequencing facility's i5/i7 index demux does *not* do — separating samples
tagged by barcodes carried on your custom PCR primers.

## Features

- **Two interfaces, one engine** — a scriptable **CLI** for pipelines and a
  modern **desktop GUI** for bench scientists (zero extra dependencies).
- **Fully configurable panel** — anchor sequences, barcode length, barcode
  lists, and mismatch policy are all user-defined; nothing is hard-coded.
- **Robust matching** — Hamming (default) or Levenshtein matching, configurable
  mismatch budget, and paired-end reconciliation that refuses to guess.
- **Transparent failures** — unassigned reads are counted by reason
  (`no_anchor`, `barcode_unmatched`, `barcode_ambiguous`, `paired_end_disagreement`).
- **Flexible I/O** — FASTQ and/or FASTA output, gzip-aware, optional per-strand
  (`_fwd`/`_rev`) splitting, streaming writes for bounded memory.
- **Fast** — optional `isal` acceleration; validated on real 6.6 M-read libraries.

## Install

```bash
pip install ampliscan[fast,excel]     # 'fast' = isal (2–5× faster gzip)
```

Or from source:

```bash
git clone https://github.com/Ayushmania2002/ampliscan
cd ampliscan
pip install -e .[fast,excel]
```

## Usage

### Desktop GUI (no terminal needed)

```bash
ampliscan gui          # or just: ampliscan
```

Configure the panel (or load it from YAML/Excel), pick your R1/R2 files and an
output folder, set options, and run — with a live progress bar, a reads-per-bin
heatmap, and an "open output folder" button.

> On Windows, `ampliscan-gui` (and double-clicking it) opens with no console
> window. On Debian/Ubuntu, tkinter may need `sudo apt install python3-tk`.

### Command line

```bash
ampliscan panel  Data.xlsx --out panel.yaml          # Excel -> panel YAML
ampliscan validate panel.yaml                        # lint it
ampliscan sniff  --r1 R1.fq.gz --r2 R2.fq.gz --config panel.yaml   # quick preview
ampliscan demux  --r1 R1.fq.gz --r2 R2.fq.gz --config panel.yaml --out bins/ --format both
```

Useful `demux` flags: `--format {fastq,fasta,both}`, `--split-by-strand`
(separate `_fwd`/`_rev` files), `--workers N` (CPU cores; see the speed note).

### Python API

```python
from ampliscan import Panel, MatchPolicy, demux_to_bins

panel = Panel.from_yaml("panel.yaml")
stats = demux_to_bins("R1.fq.gz", panel, out_dir="bins/",
                      r2_path="R2.fq.gz",
                      policy=MatchPolicy(barcode_mismatch=1),
                      output_format="both")
print(f"{stats.assigned_fraction:.1%} assigned")
```

See [`examples/ampliscan_demo.ipynb`](examples/ampliscan_demo.ipynb) for an
interactive walkthrough.

## A note on speed

For gzipped Illumina input the job is **I/O-bound**, not CPU-bound: gzip
(de)compression dominates and the barcode matching is cheap.

- Installing **`isal`** (`pip install ampliscan[fast]`) is the real speedup
  (~1.7× on a full pool) — used automatically when present.
- `--workers > 1` does **not** help gzipped input; raise it only for
  uncompressed input, `--indels` matching, or very large panels.

## How it works

See the accompanying paper, [`PAPER.md`](PAPER.md), for the method, the
octanucleotide rationale, validation, and benchmarks.

## Citing ampliscan

If you use ampliscan in your research, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Ayushman Mallick (2026). *ampliscan: an anchor-based amplicon demultiplexer
> with a desktop GUI for barcode-tagged sequencing reads* (Version 0.1.0)
> [Software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## License

[MIT](LICENSE) © 2026 Ayushman Mallick.
