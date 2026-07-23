<!--
  ┌─────────────────────────────────────────────────────────────────────┐
  │  BEFORE YOU UPLOAD TO ZENODO, FILL IN THE PLACEHOLDERS MARKED [...]  │
  │  - Author full name, affiliation, ORCID                             │
  │  - Zenodo DOI (Zenodo mints it automatically when you publish;       │
  │    paste it back here + into CITATION once you have it)             │
  │  - GitHub URL once the repo is created                              │
  └─────────────────────────────────────────────────────────────────────┘
-->

# ampliscan: an anchor-based amplicon demultiplexer with a desktop GUI for barcode-tagged sequencing reads

**Author:** Ayushman Mallick<sup>1</sup>

<sup>1</sup> [Department / Laboratory], [Institution], [City, Country]

**Correspondence:** ayushmania2002@gmail.com · **ORCID:** [0000-0000-0000-0000]

**Version:** 0.1.0 · **License:** MIT · **DOI:** [10.5281/zenodo.XXXXXXX] · **Repository:** [https://github.com/<user>/ampliscan]

---

## Abstract

High-throughput amplicon sequencing routinely pools many samples into a single
run and distinguishes them afterwards using short synthetic barcodes.
When those barcodes are carried on custom primers and sit *inside* the
sequenced fragment — immediately adjacent to a constant "anchor" sequence —
standard index-based demultiplexing performed by the sequencing facility does
not separate them. **ampliscan** is a lightweight, dependency-minimal Python
package that performs this second, inner level of demultiplexing. It locates a
user-defined pair of constant anchors, reads the variable barcode adjacent to
each, and assigns every read (or read pair) to the corresponding sample bin.
By matching the anchor and barcode together as a single **octanucleotide**
unit, ampliscan tolerates sequencing errors within the barcode while keeping
misassignment low. The tool exposes the same engine through a scriptable
command-line interface and a modern, zero-dependency desktop GUI, so it is
usable by both pipeline developers and bench scientists. On a synthetic
benchmark ampliscan recovers 100 % of error-free reads to their correct sample
and ~98 % under a realistic 0.5 % per-base error rate; on a real 6.6-million-read
paired-end library it assigned 85.8 % of reads. ampliscan is released under the
MIT license and archived on Zenodo for citation.

**Keywords:** amplicon sequencing; demultiplexing; barcodes; primers; Illumina;
bioinformatics; FASTQ; Python; GUI

---

## 1. Introduction

Sequencing many biological samples on one flow cell is far cheaper than running
them individually. To make this possible, each sample is tagged with a unique
DNA **barcode** before the samples are pooled ("multiplexed") and sequenced
together. After sequencing, the pooled reads must be sorted back into
per-sample groups — a step called **demultiplexing** — using the barcode on
each read as a label [1].

Illumina instruments and their associated software perform a first level of
demultiplexing automatically, using the i5/i7 **index** reads that sit in the
sequencing adapters. However, many amplicon experiments add a *second*,
in-house layer of barcoding: short oligonucleotides built into the forward and
reverse PCR primers, positioned inside the amplified fragment. Because these
inner barcodes are not index reads, the facility's pipeline leaves them intact,
and the returned FASTQ files still contain a mixture of samples. Separating
them is the responsibility of the researcher.

This inner barcode is typically flanked by a short **constant sequence** shared
by every primer — an *anchor*. A representative fragment has the structure:

```
[P5 adapter] · ANCHOR5 · [F-barcode] · <target insert> · [R-barcode] · ANCHOR3 · [P7 adapter]
             └──── octanucleotide ────┘                  └──── octanucleotide ────┘
```

ampliscan is a small, focused tool for exactly this task: demultiplexing reads
by inner, anchor-flanked barcodes, with an emphasis on configurability,
correctness, honest reporting of failures, and accessibility to users who do
not work at the command line.

## 2. Statement of need

General adapter- and barcode-handling tools such as cutadapt are powerful and
widely used [2], but they are command-line-only and their flexibility can make
the specific "anchor + inner barcode" workflow verbose to express, especially
for users without scripting experience. In practice, many wet-lab researchers
receive their pooled FASTQ files and have no straightforward, point-and-click
way to split them by their custom primer barcodes.

ampliscan addresses this gap with three design goals:

1. **A single, explicit model of the barcoding scheme.** The user declares the
   two anchors, the barcode length, and the forward/reverse barcode tables in a
   small human-readable configuration (or interactively). Nothing is hard-coded.
2. **Two front-ends over one engine.** The identical demultiplexing core is
   reachable from a scriptable command-line interface (for pipelines) and from
   a desktop graphical interface (for bench scientists). The GUI requires no
   dependencies beyond the Python standard library, so `pip install ampliscan`
   is sufficient.
3. **Transparency about what could not be assigned.** Every read that is not
   placed into a sample bin is counted under an explicit reason, so failures in
   library preparation or sequencing are visible rather than silently dropped.

## 3. The octanucleotide barcoding scheme

The central methodological choice in ampliscan is to treat the constant anchor
and the variable barcode as one unit. In the common design where both are four
nucleotides long, this unit is an **octanucleotide** (e.g. `gcttGCGT`).

Matching on the barcode alone is fragile: a single sequencing error in a 4-nt
barcode can either push a read into the wrong sample's bin (a silent error) or
cause it to be discarded (lost data). The constant anchor mitigates both
problems. Because the anchor is identical in every read, it acts as a fixed
landmark that pins the alignment to the correct position; the variable barcode
immediately adjacent to it is then read off and matched. This is the same
principle that underlies anchored adapter trimming in established tools [2], and
it materially improves discrimination in our benchmarks (Section 5).

## 4. Implementation

ampliscan is written in Python (≥ 3.9). The core engine depends only on the
standard library plus PyYAML; optional extras enable Excel panel parsing
(pandas/openpyxl), accelerated compression (`isal`), and richer notebook
visualisation. The package is organised into small, single-responsibility
modules: configuration and validation, anchor localisation, approximate barcode
matching, the per-read demultiplexing pipeline, streaming input/output, a
synthetic-read generator for testing, an Excel-panel importer, a parallel
file-to-bins driver, a command-line interface, and a Tkinter GUI.

### 4.1 Panel definition and validation

A *panel* specifies the two forward-strand anchors, the barcode length, and the
forward and reverse barcode dictionaries; the reverse-strand anchors default to
the reverse complement of the forward anchors but may be overridden for
asymmetric designs. On construction the panel is validated: barcode lengths
must be consistent, sequences must be valid DNA, and barcodes within each set
must be unique. Panels can be written to and read from YAML, or generated
automatically from a spreadsheet of example amplicons.

### 4.2 Anchor localisation

For each read the 5′ anchor is located by a left-to-right search and the 3′
anchor by a **right-to-left** search. Searching the 3′ anchor from the end is a
deliberate and important detail: because anchors are short, a matching sequence
can occur by chance inside the biological insert, and a naive left-to-right
search would frequently lock onto that spurious internal copy rather than the
true flanking anchor. In development, correcting this single behaviour raised
the fraction of correctly assigned synthetic reads from ~81 % to ~100 %.
Anchor matching tolerates a configurable number of substitutions (default 0).

### 4.3 Approximate barcode matching

Once the anchors are found, the flanking barcodes are extracted and matched
against the panel. Matching uses Hamming distance by default (substitutions
only), with an optional banded Levenshtein mode (insertions and deletions) for
indel-prone chemistries. A read is assigned only when exactly one panel barcode
lies within the mismatch budget (default 1); a barcode that is equidistant from
two panel entries is reported as **ambiguous** rather than guessed.

### 4.4 Paired-end reconciliation and strand handling

For paired-end input, R1 and R2 are demultiplexed independently and their
sample assignments are then required to agree. When they disagree — indicating
uncorrectable noise in one mate — the read pair is withheld rather than assigned
to a possibly-wrong sample, on the principle that silent cross-sample
contamination is worse than a modest loss of data. Each read also records the
orientation in which it matched: the `+` strand (forward anchors, R1-like) or
the `−` strand (reverse-complement anchors, R2-like). Reads may optionally be
written to separate `_fwd`/`_rev` files per sample.

### 4.5 Classification of unassigned reads

Reads that cannot be placed are counted under one of four explicit reasons:
*no_anchor* (neither flanking anchor found), *barcode_unmatched* (barcode too
distant from any panel entry), *barcode_ambiguous* (equidistant from two
entries), and *paired_end_disagreement* (mates assigned to different samples).
This breakdown makes systematic wet-lab or sequencing problems visible.

### 4.6 Input, output, and performance

ampliscan reads FASTQ and FASTA, plain or gzip-compressed, and auto-detects the
format from the file extension. Output is written incrementally through a
streaming per-bin writer, so peak memory stays bounded (approximately one work
chunk) regardless of total input size; one file is produced per sample, and, if
requested, in more than one format simultaneously (FASTQ and FASTA).

A notable and deliberately reported empirical finding is that, for typical
gzip-compressed Illumina input, the workload is **I/O-bound rather than
CPU-bound**: decompression and compression dominate, while the barcode matching
itself is inexpensive. Consequently, distributing the matching across multiple
processes yields little benefit and can even be slower, because serialising
millions of records between processes costs more than the matching it
parallelises. The effective optimisation is faster (de)compression: installing
the optional `isal` library reduces the runtime of a full 6.6-million-read
paired-end demultiplexing from ~181 s to ~104 s on a single core (Section 5.3).
A user-selectable worker count is nonetheless provided for the cases where
matching *is* the bottleneck (uncompressed input, Levenshtein matching, or very
large barcode panels), with single-process operation as the default.

### 4.7 Interfaces

The command-line interface provides subcommands to build a panel from a
spreadsheet, validate a panel, preview the bin distribution on a sample of
reads, and run a full demultiplexing. The desktop GUI, built with the standard
library's Tkinter toolkit, lets the user enter or load a panel, choose input
and output locations, set matching options, choose output formats, select a
worker count, and run — displaying a live progress indicator, per-sample read
counts as a heatmap, and the breakdown of unassigned reads. Because the GUI
carries no third-party dependencies, it is available immediately after
installation and can be packaged as a standalone executable for users without a
Python environment.

## 5. Validation and results

### 5.1 Synthetic benchmark

To validate correctness under controlled conditions, ampliscan includes a
generator that synthesises paired-end reads for every barcode combination in a
panel, embedding an adjustable per-base substitution error and encoding each
read's true sample of origin in its identifier. Assignment accuracy can then be
scored exactly.

With no sequencing error, ampliscan assigns 100 % of reads to their correct
sample. Under a realistic 0.5 % per-base error rate and a mismatch budget of 1,
it assigns ~98 % of reads, with the small remainder withheld for legitimate
reasons (uncorrectable barcodes, ambiguity, or paired-end disagreement) rather
than misassigned. These results are encoded as automated tests that accompany
the package.

### 5.2 Real dataset

On a real paired-end amplicon library of 6,592,457 read pairs from an Illumina
NovaSeq X run — a pool for which the facility had already performed i5/i7 index
demultiplexing — ampliscan assigned 85.8 % of reads to the 20 expected
forward × reverse sample bins using a four-nucleotide panel with anchors
`GCTT`/`ACAG`. The remaining reads were accounted for transparently: ~8.2 %
lacked a detectable anchor, ~5.5 % carried an unmatched barcode, and small
fractions were ambiguous or showed paired-end disagreement. This breakdown is
itself a useful quality-control readout for the library preparation.

### 5.3 Runtime

On the same 6.6-million-read dataset, single-core throughput was ~36,000
read pairs · s⁻¹ using the standard-library gzip implementation (≈ 181 s total),
rising to ~63,000 read pairs · s⁻¹ with the optional `isal` accelerator (≈ 104 s
total). Increasing the worker count did not improve — and in some
configurations slightly worsened — wall-clock time, consistent with the
I/O-bound characterisation in Section 4.6. Peak memory remained modest because
output is streamed rather than accumulated.

## 6. Usage example

Command line:

```bash
pip install ampliscan[fast,excel]

# build a panel from a spreadsheet of example amplicons
ampliscan panel Data.xlsx --out panel.yaml

# demultiplex a paired-end library into per-sample FASTQ files
ampliscan demux --r1 R1.fq.gz --r2 R2.fq.gz \
                --config panel.yaml --out bins/ \
                --format both --barcode-mismatch 1
```

Python API:

```python
from ampliscan import Panel, MatchPolicy, demux_to_bins

panel = Panel.from_yaml("panel.yaml")
stats = demux_to_bins(
    "R1.fq.gz", panel, out_dir="bins/",
    r2_path="R2.fq.gz",
    policy=MatchPolicy(barcode_mismatch=1),
    output_format="both",
)
print(stats.assigned_fraction)
```

Graphical interface:

```bash
ampliscan gui      # or simply: ampliscan
```

## 7. Relationship to existing tools

ampliscan does not aim to replace comprehensive read-processing suites such as
cutadapt [2], which offer a broader range of trimming and filtering operations.
Its contribution is narrower and complementary: an explicit, validated model of
the anchor-plus-inner-barcode demultiplexing task; transparent, per-reason
accounting of unassigned reads; and — uniquely among comparable tools — a
zero-dependency graphical interface that makes the operation accessible to
researchers who do not use the command line.

## 8. Limitations and future work

ampliscan currently targets substitution-dominated short-read chemistries;
although an indel-aware matching mode is provided, it has not been extensively
validated on long-read data. It does not perform quality trimming or adapter
removal beyond the demultiplexing step, and it does not yet read or write
aligned (BAM/SAM) formats. Planned additions include an exportable HTML/PDF
quality-control report, optional per-sample consensus generation, and
alignment-format input. Contributions are welcome through the project
repository.

## 9. Availability and requirements

- **Project name:** ampliscan
- **Repository:** [https://github.com/<user>/ampliscan]
- **Archived version (this release):** DOI [10.5281/zenodo.XXXXXXX]
- **Operating systems:** platform-independent (Windows, macOS, Linux)
- **Programming language:** Python ≥ 3.9
- **Dependencies:** PyYAML (core); optional: pandas + openpyxl (Excel),
  `isal` (fast gzip), Jupyter/ipywidgets/plotly (notebook)
- **License:** MIT

## 10. How to cite

If you use ampliscan in your research, please cite this archived release:

> Ayushman Mallick (2026). *ampliscan: an anchor-based amplicon
> demultiplexer with a desktop GUI for barcode-tagged sequencing reads*
> (Version 0.1.0) [Software]. Zenodo. https://doi.org/[10.5281/zenodo.XXXXXXX]

BibTeX:

```bibtex
@software{ampliscan_2026,
  author  = {Mallick, Ayushman},
  title   = {ampliscan: an anchor-based amplicon demultiplexer with a
             desktop GUI for barcode-tagged sequencing reads},
  year    = {2026},
  version = {0.1.0},
  publisher = {Zenodo},
  doi     = {10.5281/zenodo.XXXXXXX},
  url     = {https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

## 11. Acknowledgements

The author thanks [PI / collaborators / funding source — optional] for the
sequencing data and guidance that motivated this tool.

## References

1. Illumina, Inc. *Indexed sequencing and sample demultiplexing overview.*
   Technical documentation.
2. Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput
   sequencing reads. *EMBnet.journal*, 17(1), 10–12.
   https://doi.org/10.14806/ej.17.1.200
3. Cock, P. J. A., Fields, C. J., Goto, N., Heuer, M. L., & Rice, P. M. (2010).
   The Sanger FASTQ file format for sequences with quality scores, and the
   Solexa/Illumina FASTQ variants. *Nucleic Acids Research*, 38(6), 1767–1771.
   https://doi.org/10.1093/nar/gkp1137
4. python-isal: ISA-L accelerated compression for Python.
   https://github.com/pycompression/python-isal
