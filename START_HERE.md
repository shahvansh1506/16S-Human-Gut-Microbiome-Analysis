> **GitHub export notice:** This repository is a filtered copy of the original team ZIP. It **does not include** the three trimmed pilot FASTQ files, original paper PDFs, or the 1,084 raw trimming logs. See the root `README.md` for what was retained. The original team guide below describes the complete team archive. The `pilot` mode will not work in this export until the original pilot files are restored locally.

# Final team copy · 17 November 2026

On Windows, right-click `Gut16S-Team-Final.zip` and choose **Extract All**. Open the extracted `Gut16S-Team-Final` folder. Everyone should open the project from this copy on their PC.

Open `gut16s-team.Rproj` in RStudio. Paste the following into the **Console** at the bottom left, one line at a time, and press Enter after each line.

```r
source("workflow.R")
run_team(mode = "inventory")
```

The inventory should show six rows marked `complete` and the same counts on all three PCs. It checks the saved result files and writes `results/study_inventory.tsv` for Excel. The `settings_status` column in that file says whether analysis settings were reviewed. `status` describes the result files that are present in this project. This check uses base R and needs no downloads or package installations.

### Open the saved results

```r
shell.exec(normalizePath("results/combined", mustWork = TRUE))
```

The folder contains the six-study graph, results report, tables and the one-page `Gut_Microbiome_Results_Handout.pdf`. The verified full outputs for each study are under `results/full`. The final copy contains all six studies and 1,167 processed runs. PRJEB20022 uses the forward FASTQ from each of 859 selected runs. The archive manifest still records the paired files that were deposited.

To regenerate the combined graph and results report from those saved tables, use `source("make_final_report.R")`. It requires ggplot2. The saved graph and report can be opened directly on every PC.

## Optional analysis and pilot demonstration

Each teammate can request their assigned studies in a single command. Start with `preview`, which downloads and checks up to three candidate runs per study. Choose the line with your own name.

```r
run_team("Vansh", "preview")
run_team("Bhavya", "preview")
run_team("Krupa", "preview")
```

Open the PDFs under `results/preview/<project>/`. They show per-base read quality before trimming. Where candidate selection is pending, the workflow prints that study's ID and moves to the next study. Files that match the archive's checksum are reused.

The compressed FASTQs selected for the final analysis total about 4.5 GB across all six studies. The **original team ZIP**, not this GitHub export, includes the three trimmed pilot inputs. Full input downloads and temporary processing folders stay with the original project copies.

## Software

Use R 4.5.x, Bioconductor 3.22, DADA2 1.38.0 and ggplot2. If DADA2 and ggplot2 are already installed, check them in the RStudio Console with `packageVersion("dada2")` and `packageVersion("ggplot2")`. Install missing packages with the commands below in a fresh R 4.5.x session.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("dada2", version = "3.22", update = FALSE, ask = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
```

For studies that need primer trimming, install Cutadapt 4.7. The workflow detects it on the system PATH, in a Python installation, or in an Ubuntu WSL distribution on Windows. After `source("workflow.R")`, run `source("fix_cutadapt_windows.R")` to detect Windows distributions with names such as `Ubuntu-24.04`. Teammates can follow the [Cutadapt installation documentation](https://cutadapt.readthedocs.io/en/v4.7/installation.html). If Ubuntu WSL already has Python 3 and `venv`, these commands install it in an isolated environment:

```bash
python3 -m venv ~/gut16s-cutadapt
~/gut16s-cutadapt/bin/python -m pip install cutadapt==4.7
~/gut16s-cutadapt/bin/cutadapt --version
```

Use the Ubuntu terminal for those three commands. For that layout the workflow detects `~/gut16s-cutadapt/bin/cutadapt` automatically on Windows. The [Bioconductor DADA2 1.38.0 page](https://bioconductor.org/packages/3.22/bioc/html/dada2.html) documents the R version. The reference downloads once from [Zenodo](https://zenodo.org/records/14169026) and is verified by MD5. If you already have that exact file, put it in the team's project `reference` folder before running.

## Test the workflow on a few samples

`pilot` runs the original three PRJNA687346 trimmed samples. Use it to check the R installation and the shared pipeline.

```r
run_team(mode = "pilot")
```

`trial` runs the first, middle and last eligible runs of each selected study. The same code and reference are used for the full run. Trial files are written under `results/trial/<project>/`.

```r
run_team("PRJNA668251", "trial")
```

To explain how settings were chosen, open a `read_tracking.tsv` file. Compare `filtered` to `input`, `merged` to `inferred` for paired reads, and `bacterial` to `nonchim`. A sample with zero bacterial reads stays in the tracking table and has no composition percentages. `input` counts reads or read pairs entering DADA2 after Cutadapt, when applicable. Files under `trimming_logs` show how many reads matched the primers.

The paired trimming configuration keeps read pairs with both primer matches. The three-sample knee OA pilot uses its original trimming rule, which kept reads even if Cutadapt did not find the primer. The audit records where the archives appear to have removed primers or combined the two reads before deposit.

## Run the full studies and combine the results

The workflow uses `runs.tsv` for specimen IDs, candidate selection, archive layout and sequencing batches. The final settings in `studies.tsv` match the verified runs. `status = ready` in that file means the settings and selection were reviewed. `include = yes` selects a run. `trim_mode = forward` uses the forward file only. The audit records the six study decisions. Pilot and trial results use their own folders.

Full mode first checks the saved results. If a study needs processing, `max_download_gb` limits the sum of its selected compressed FASTQ sizes. Leave additional free space for trimmed reads, filtered reads and results. Here are the commands for each student:

```r
run_team("Vansh", "full", max_download_gb = 3)
run_team("Bhavya", "full", max_download_gb = 3)
run_team("Krupa", "full", max_download_gb = 8)
```

The combined graph, study means, and sample percentages are in `results/combined/`. Each sample contributes equally to its own study's mean. Per-study ASVs and taxonomy are available as RDS files. Interpret the plot as descriptive because study cohorts and laboratory protocols differ.

Results from the shared ZIP can be reused when the archived checksums, recorded inputs and settings match. For a fresh run, use a separate copy of the project and keep the verified results for comparison. A failed run prints the working folder path. Raw downloads with valid checksums will be reused on retry.
