# Validation

- The packaged project loaded in R 4.5.2 and its inventory listed all six projects and 1,167 candidate runs.
- The three original trimmed PRJNA687346 samples ran through DADA2 and SILVA taxonomy. Nonchimera read totals matched the original analysis at 39,158, 38,125 and 32,724. Both analyses yielded 448 ASVs.
- Downloaded and verified the ENA checksums for both mates of three PRJNA668251 runs. Cutadapt 4.7 ran through Ubuntu WSL and retained 98.9–99.2% of pairs with both primer matches.
- The quality preview generated a PDF for those three pairs using the Windows R installation.
- A test on 10,000 paired reads from PRJNA668251 retained 9,338 pairs through filtering and 8,356 through denoising and merging. A separate filtering check on the entire same sample retained 122,232 of 129,939 pairs entering DADA2. The paired tests used the first ENA read-header field to identify mates.
- The combined table and plot passed a two-study test with known percentages and a phylum absent in one study. Each study mean weighted its retained samples equally.

The final team copy includes verified full results for six projects and 1,167 eligible runs. The run IDs match the archive selection, the phylum percentages match the study results, and the settings match the original completion records. The release fingerprints in `verified_full_results.tsv` tie the files in the ZIP to these verified copies. PRJEB20022 contains the 859 selected archived subject runs processed with forward reads only. See `study_audit.md` for the scope of each study and the limits on cross-study interpretation.
