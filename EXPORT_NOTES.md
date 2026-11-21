# Repository export notes

This folder was prepared from the uploaded `Gut16S-Team-Final.zip` and `Six_Selected_Papers.zip` files.

## Preserved

- The original R scripts, `.Rproj` metadata, study and run inventory tables.
- Six saved study output sets: read tracking, per-sample phylum summaries, ASV tables, taxonomy assignments, run selections, and software session records.
- Team-generated combined plots, reports, examples, audit notes, and reported validation checks.
- All original 42 checksum-matched files listed in `verified_full_results.tsv` were retained without modification.

## Intentionally excluded

- Three trimmed public-source FASTQ inputs (human-associated microbial sequence reads).
- 1,084 `cutadapt` trimming logs; the summary tables and reproducibility metadata remain.
- Six third-party full-text PDF research papers. The README/audit link to the primary publications instead.
- Any large reference database downloads and temporary run intermediates.

## Validation performed on this export

- Listed files in both supplied ZIPs and inspected the scripts, result tables, original team reports, validation notes, and the supporting paper titles.
- Matched 42 per-study file MD5 hashes against the supplied `verified_full_results.tsv` manifest.
- Confirmed that archived `results/combined/study_read_summary.tsv` reports 1,167 runs across six studies.
- Did not independently install or run R/DADA2 or download/reprocess raw FASTQ data in this environment.

## Before public release

1. Check with Bhavya and Krupa that they agree to publish the jointly created code and archived figures/tables with appropriate credit.
2. Choose a repository license only after the team agrees.
3. Optionally update the contributor display names to full names if everyone wants that.
