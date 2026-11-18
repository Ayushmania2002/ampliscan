# Changelog

## 1.2.2
- Fix: the 3' anchor search is now bounded at the Illumina adapter (`MatchPolicy.adapter_trim`, default `AGATCGGAAGAGC`; CLI `--adapter SEQ` / `--no-adapter-trim`). Reads that run through the insert into the adapter and sample index could match a stray copy of the anchor inside the index (e.g. index `AGTATCTTGC` gives `...CCAGTCACAGTATC...`), so the R barcode was read from the wrong place and the read was lost as `barcode_unmatched`.
- Fix: when neither mate of a pair is assigned, the reported reason is now the more informative one instead of always R2's. A failed (all-N) Read 2 no longer masks R1's real reason as `no_anchor`.

## 1.2.1
- Reads with no insert between the barcodes are no longer written to sample bins (new unassigned reason `target_too_short`, `--min-target-length`).

## 1.2.0
- `ampliscan demux --report` writes a self-contained HTML report.

## 1.1.0
- GUI shows two heatmaps, R1 alone vs R2 alone.
