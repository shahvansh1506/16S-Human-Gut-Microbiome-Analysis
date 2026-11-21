# Shared human gut 16S workflow. Open the project before sourcing this file.

reference_file <- file.path("reference", "silva_nr99_v138.2_toGenus_trainset.fa.gz")
reference_url <- paste0("https://zenodo.org/records/14169026/files/", basename(reference_file), "?download=1")
reference_md5 <- "1764e2a36b4500ccb1c7d5261948a414"
pilot_ids <- c("SRR13298477", "SRR13298480", "SRR13298496")

read_table <- function(path) {
  read.delim(path, quote = "", comment.char = "", colClasses = "character",
             na.strings = character(), check.names = FALSE)
}

write_table <- function(data, path) {
  write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

fetch_file <- function(url, path, md5, bytes = NA_real_) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  valid <- function(file) {
    file.exists(file) && file.size(file) > 0 &&
      (is.na(bytes) || file.size(file) == bytes) &&
      identical(tolower(unname(tools::md5sum(file))), tolower(md5))
  }
  if (valid(path)) return(path)
  tmp <- tempfile("download-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  old <- options(timeout = max(180, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  for (attempt in seq_len(3)) {
    message("Downloading ", basename(path), " (attempt ", attempt, ")")
    ok <- tryCatch(download.file(url, tmp, mode = "wb", method = "libcurl", quiet = TRUE) == 0L,
                   error = function(e) { message(conditionMessage(e)); FALSE })
    if (ok && valid(tmp)) {
      if (file.exists(path)) stop("Move the invalid cached file before retrying: ", path)
      if (!file.rename(tmp, path)) stop("Could not save ", path)
      return(path)
    }
  }
  stop("Download or checksum failed: ", url)
}

load_inputs <- function() {
  studies <- read_table("studies.tsv")
  runs <- read_table("runs.tsv")
  stopifnot(!anyDuplicated(studies$project), !anyDuplicated(runs$run),
            all(runs$project %in% studies$project),
            all(grepl("^PRJ[A-Z]+[0-9]+$", studies$project)),
            all(grepl("^[SED]RR[0-9]+$", runs$run)))
  list(studies = studies, runs = runs)
}

runs_for_study <- function(study, runs) {
  runs <- runs[runs$project == study$project & runs$include == "yes", , drop = FALSE]
  if (nrow(runs) && study$trim_mode == "forward") {
    runs$layout <- rep("SINGLE", nrow(runs))
    runs[c("url_r", "md5_r", "bytes_r")] <- ""
  }
  runs[order(runs$run), , drop = FALSE]
}

full_result_files <- c("complete.rds", "used_runs.tsv", "read_tracking.tsv", "phylum_percent.tsv",
                       "asv_counts.rds", "taxonomy.rds", "session_info.txt")

full_result_available <- function(study, runs) {
  folder <- file.path("results", "full", study$project)
  paths <- file.path(folder, full_result_files)
  if (!nrow(runs) || !all(file.exists(paths) & file.size(paths) > 0)) return(FALSE)
  record <- tryCatch(readRDS(paths[1]), error = function(e) NULL)
  if (!is.list(record) || !identical(record$scope, "full") || !identical(record$ref, reference_md5) ||
      !identical(record$dada2, "1.38.0")) return(FALSE)

  fields <- c("project", "region", "primer_f", "primer_r", "trim_mode",
              "trunc_f", "trunc_r", "maxee_f", "maxee_r")
  settings_match <- all(vapply(fields, function(field) identical(record$study[[field]], study[[field]]), logical(1)))
  error_model <- function(x) if (length(x) && nzchar(x)) x else "default"
  if (!settings_match || !identical(error_model(record$study$error_model), error_model(study$error_model))) return(FALSE)
  if (!identical(names(runs), names(record$runs)) || !setequal(runs$run, record$runs$run)) return(FALSE)
  original_order <- match(runs$run, record$runs$run)
  if (!all(vapply(names(runs), function(key) identical(runs[[key]], record$runs[[key]][original_order]), logical(1)))) return(FALSE)

  # Preserve the original result records when the team distributes verified runs.
  # Other results must have been produced by this exact workflow script.
  archive <- "verified_full_results.tsv"
  if (file.exists(archive)) {
    published <- read_table(archive)
    published <- published[published$project == study$project, , drop = FALSE]
    if (nrow(published)) {
      if (nrow(published) != length(full_result_files) ||
          !setequal(published$file, full_result_files)) return(FALSE)
      fingerprints <- unname(tools::md5sum(paths))
      if (identical(published$md5[match(full_result_files, published$file)], fingerprints)) return(TRUE)
    }
  }
  identical(record$workflow, unname(tools::md5sum("workflow.R")))
}

inventory <- function() {
  data <- load_inputs()
  d <- data$studies[c("owner", "project", "region", "status", "note")]
  d$settings_status <- d$status
  d$runs_in_archive <- vapply(d$project, function(p) sum(data$runs$project == p), integer(1))
  d$candidate_runs <- vapply(d$project, function(p) sum(data$runs$project == p & data$runs$include == "yes"), integer(1))
  selected <- lapply(seq_len(nrow(d)), function(i) runs_for_study(data$studies[i, , drop = FALSE], data$runs))
  complete <- vapply(seq_len(nrow(d)), function(i) full_result_available(data$studies[i, , drop = FALSE], selected[[i]]), logical(1))
  folders <- dir.exists(file.path("results", "full", d$project))
  d$status[folders & !complete] <- "check_files"
  d$status[complete] <- "complete"
  d$completed_runs <- ifelse(complete, d$candidate_runs, 0L)
  d$fastq_GB <- vapply(selected, function(rows) {
    sum(as.numeric(c(rows$bytes_f, rows$bytes_r)), na.rm = TRUE) / 1e9
  }, numeric(1))
  d$review_runs <- vapply(d$project, function(p) sum(data$runs$project == p & data$runs$include == "review"), integer(1))
  print(d[c("owner", "project", "status", "runs_in_archive", "candidate_runs", "completed_runs")], row.names = FALSE)
  dir.create("results", showWarnings = FALSE)
  write_table(d, file.path("results", "study_inventory.tsv"))
  invisible(d)
}

download_reads <- function(runs) {
  out <- data.frame(f = character(nrow(runs)), r = character(nrow(runs)))
  for (i in seq_len(nrow(runs))) {
    row <- runs[i, ]
    for (mate in if (row$layout == "PAIRED") c("f", "r") else "f") {
      url <- row[[paste0("url_", mate)]]
      md5 <- row[[paste0("md5_", mate)]]
      if (!grepl("^https://ftp\\.sra\\.ebi\\.ac\\.uk/vol1/fastq/", url) ||
          !grepl("^[a-fA-F0-9]{32}$", md5)) stop("No verified ENA file route for ", row$run)
      out[i, mate] <- fetch_file(url, file.path("data", "raw", row$project,
                                paste0(row$run, "_", mate, ".fastq.gz")),
                                md5, as.numeric(row[[paste0("bytes_", mate)]]))
    }
  }
  out
}

cutadapt_runner <- function() {
  candidates <- list(c(Sys.which("cutadapt")), c(Sys.which("python"), "-m", "cutadapt"),
                     c(Sys.which("py"), "-m", "cutadapt"))
  if (.Platform$OS.type == "windows" && nzchar(Sys.which("wsl.exe"))) {
    wsl <- Sys.which("wsl.exe")
    wsl_home <- tryCatch(suppressWarnings(system2(wsl, c("-d", "Ubuntu", "--", "printenv", "HOME"),
                                                 stdout = TRUE, stderr = TRUE)),
                         error = function(e) character())
    wsl_prefix <- c(wsl, "-d", "Ubuntu", "--cd", getwd(), "--")
    candidates <- c(candidates, list(c(wsl_prefix, "cutadapt")))
    if (length(wsl_home) == 1L && grepl("^/[^\r\n]+$", wsl_home)) {
      candidates <- c(candidates, list(c(wsl_prefix, paste0(wsl_home, "/gut16s-cutadapt/bin/cutadapt"))))
    }
  }
  for (candidate in candidates) {
    if (!nzchar(candidate[1])) next
    runner <- local({
      command <- candidate[1]
      prefix <- candidate[-1]
      function(args) {
        if (basename(command) == "wsl.exe") args <- chartr("\\", "/", args)
        arguments <- c(prefix, args)
        quote <- !grepl("^[A-Za-z0-9_./:=+-]+$", arguments)
        arguments[quote] <- shQuote(arguments[quote])
        suppressWarnings(system2(command, arguments, stdout = TRUE, stderr = TRUE))
      }
    })
    version <- tryCatch(runner("--version"), error = function(e) character())
    if (is.null(attr(version, "status")) &&
        any(trimws(version) == "4.7")) return(runner)
  }
  stop("Install Cutadapt 4.7 or make it available in Ubuntu WSL. See START_HERE.md.")
}

preview_project <- function(runs, folder) {
  reads <- download_reads(runs)
  # Quality previews sample only 10,000 reads per file. Avoid starting a pool of
  # R workers on each file, which is costly on Windows laptops.
  previous_workers <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(previous_workers), add = TRUE)
  grDevices::pdf(file.path(folder, "read_quality.pdf"), width = 11, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(dada2::plotQualityProfile(reads$f, n = 10000))
  if (any(nzchar(reads$r))) print(dada2::plotQualityProfile(reads$r[nzchar(reads$r)], n = 10000))
  write_table(runs[c("project", "run", "sample", "layout", "batch")], file.path(folder, "preview_runs.tsv"))
  message("Quality plots saved in ", folder)
}

trim_read <- function(reads, study, run, folder, runner) {
  if (study$trim_mode == "already_trimmed") return(reads)
  if (!study$trim_mode %in% c("pilot", "both", "forward")) stop("Review primer handling for ", run)
  paired <- nzchar(reads$r)
  f <- file.path(folder, paste0(run, "_f.fastq.gz"))
  r <- if (paired) file.path(folder, paste0(run, "_r.fastq.gz")) else ""
  args <- c("-g", study$primer_f, "-O", nchar(study$primer_f), "-o", f)
  if (paired) {
    if (study$trim_mode == "both") {
      args <- c(args, "-G", study$primer_r, "--discard-untrimmed", "--pair-filter=any")
    }
    args <- c(args, "-p", r)
  }
  report <- runner(c(args, reads$f, if (paired) reads$r))
  writeLines(report, file.path(folder, paste0(run, "_cutadapt.log")))
  if (!is.null(attr(report, "status")) || !all(file.exists(c(f, if (paired) r)))) {
    stop("Cutadapt failed for ", run, ". See ", folder)
  }
  data.frame(f = f, r = r)
}

infer_batch <- function(runs, reads, study, folder) {
  paired <- identical(unique(runs$layout), "PAIRED")
  f <- setNames(file.path(folder, paste0(runs$run, "_f.fastq.gz")), runs$run)
  r <- setNames(file.path(folder, paste0(runs$run, "_r.fastq.gz")), runs$run)
  trunc <- as.integer(c(study$trunc_f, if (paired) study$trunc_r))
  ee <- as.numeric(c(study$maxee_f, if (paired) study$maxee_r))
  stopifnot(all(is.finite(trunc) & trunc >= 0), all(is.finite(ee) & ee > 0))
  args <- list(fwd = reads$f, filt = f, truncLen = trunc, maxEE = ee, maxN = 0,
               truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  # ENA puts the shared archive spot identifier in field 1. Instrument headers
  # in later fields may end in /1 and /2 and cannot be matched literally.
  if (paired) args <- c(args, list(rev = reads$r, filt.rev = r, matchIDs = TRUE, id.field = 1L))
  qc <- do.call(dada2::filterAndTrim, args)
  if (!identical(rownames(qc), basename(reads$f))) {
    # The filter output uses input filenames as row names.
    qc <- qc[match(basename(reads$f), rownames(qc)), , drop = FALSE]
  }
  stopifnot(!anyNA(qc), nrow(qc) == nrow(runs))
  track <- data.frame(run = runs$run, sample = runs$sample, batch = runs$batch,
                      input = qc[, "reads.in"], filtered = qc[, "reads.out"],
                      inferred = 0, merged = 0, nonchim = 0)
  ok <- track$filtered > 0
  if (!any(ok)) stop("No reads passed filtering for ", study$project, " / ", runs$batch[1])
  error_model <- study$error_model
  if (is.null(error_model) || !nzchar(error_model)) error_model <- "default"
  estimate <- switch(error_model,
                     default = dada2::loessErrfun,
                     constant_quality = dada2::noqualErrfun,
                     stop("Unknown error_model for ", study$project, ": ", error_model))
  if (error_model == "constant_quality") {
    message("Learning error rates without quality-score differences for ", study$project, ".")
  }
  err_f <- dada2::learnErrors(f[ok], multithread = FALSE, errorEstimationFunction = estimate)
  if (paired) err_r <- dada2::learnErrors(r[ok], multithread = FALSE, errorEstimationFunction = estimate)
  variants <- setNames(vector("list", sum(ok)), runs$run[ok])
  for (i in which(ok)) {
    derep_f <- dada2::derepFastq(f[i])
    inferred_f <- dada2::dada(derep_f, err = err_f, multithread = FALSE)
    track$inferred[i] <- sum(dada2::getUniques(inferred_f))
    if (paired) {
      derep_r <- dada2::derepFastq(r[i])
      inferred_r <- dada2::dada(derep_r, err = err_r, multithread = FALSE)
      variants[[runs$run[i]]] <- dada2::mergePairs(inferred_f, derep_f, inferred_r, derep_r)
    } else variants[[runs$run[i]]] <- inferred_f
    track$merged[i] <- sum(dada2::getUniques(variants[[runs$run[i]]]))
  }
  if (!any(track$merged > 0)) stop("No sequences after inference or paired read merging in ", runs$batch[1])
  variants <- variants[track$run[track$merged > 0]]
  list(table = dada2::makeSequenceTable(variants), track = track)
}

study_results <- function(study, runs, output, local_trimmed = NULL) {
  paired <- unique(runs$layout)
  if (length(paired) != 1 || !paired %in% c("SINGLE", "PAIRED")) stop("Mixed or unsupported layout in ", study$project)
  if (any(!nzchar(runs$batch)) || anyDuplicated(runs$sample)) stop("Review batches and duplicate sample IDs for ", study$project)
  if (is.null(local_trimmed) && study$trim_mode != "already_trimmed") {
    if (!nzchar(study$primer_f) || (paired == "PAIRED" && study$trim_mode == "both" && !nzchar(study$primer_r))) {
      stop("Primer settings missing for ", study$project)
    }
    runner <- cutadapt_runner()
  }
  final <- file.path(output, study$project)
  signature <- list(study = study, runs = runs, ref = reference_md5, scope = basename(output),
                    workflow = unname(tools::md5sum("workflow.R")),
                    dada2 = as.character(packageVersion("dada2")),
                    local_reads = if (length(local_trimmed)) tools::md5sum(local_trimmed) else NULL)
  if (file.exists(file.path(final, "complete.rds"))) {
    if (!identical(readRDS(file.path(final, "complete.rds")), signature)) {
      stop("Settings changed. Move the previous result folder before running: ", final)
    }
    message("Using completed results for ", study$project)
    return(invisible(final))
  }
  if (dir.exists(final)) stop("Move the incomplete result folder before retrying: ", final)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  stage <- tempfile(paste0(study$project, "-"), tmpdir = output)
  dir.create(stage)
  dir.create(file.path("cache", "work"), recursive = TRUE, showWarnings = FALSE)
  work <- tempfile(paste0(study$project, "-"), tmpdir = file.path("cache", "work"))
  dir.create(work)
  message("Processing ", study$project, " in ", stage)
  if (is.null(local_trimmed)) {
    reads <- download_reads(runs)
    if (study$trim_mode != "already_trimmed") {
      reads <- do.call(rbind, lapply(seq_len(nrow(runs)), function(i) {
        trim_read(reads[i, ], study, runs$run[i], work, runner)
      }))
    }
  } else reads <- data.frame(f = local_trimmed, r = "")
  filtered <- file.path(work, "filtered")
  dir.create(filtered)
  set.seed(123)
  batches <- lapply(split(seq_len(nrow(runs)), runs$batch), function(index) {
    infer_batch(runs[index, ], reads[index, ], study, filtered)
  })
  tables <- lapply(batches, `[[`, "table")
  raw_tab <- if (length(tables) == 1L) tables[[1L]] else dada2::mergeSequenceTables(tables = tables)
  counts <- dada2::removeBimeraDenovo(raw_tab, method = "consensus", multithread = FALSE)
  track <- do.call(rbind, lapply(batches, `[[`, "track"))
  track$nonchim <- as.numeric(rowSums(counts)[match(track$run, rownames(counts))])
  track$nonchim[is.na(track$nonchim)] <- 0
  write_table(track, file.path(stage, "read_tracking.tsv"))
  if (!ncol(counts)) stop("No ASVs retained for ", study$project)
  set.seed(123)
  tax <- dada2::assignTaxonomy(counts, reference_file, minBoot = 50, multithread = FALSE)
  tax <- tax[match(colnames(counts), rownames(tax)), , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(tax)))
  organelle <- apply(tax, 1, function(x) any(tolower(x) %in% c("chloroplast", "mitochondria", "mitochondrion")))
  bacteria <- tax[, "Kingdom"] %in% "Bacteria" & !organelle
  phylum <- tax[bacteria, "Phylum"]
  phylum[is.na(phylum)] <- "Unclassified phylum"
  if (!length(phylum)) stop("No bacterial ASVs retained for ", study$project)
  phylum_counts <- rowsum(t(counts[, bacteria, drop = FALSE]), group = phylum)
  valid <- colSums(phylum_counts) > 0
  if (!any(valid)) stop("No sample had retained bacterial reads")
  phylum_counts <- phylum_counts[, valid, drop = FALSE]
  pct <- as.data.frame(as.table(100 * sweep(phylum_counts, 2, colSums(phylum_counts), "/")))
  names(pct) <- c("phylum", "run", "percent")
  pct$project <- study$project
  pct$sample <- runs$sample[match(pct$run, runs$run)]
  track$bacterial <- as.numeric(colSums(phylum_counts)[match(track$run, colnames(phylum_counts))])
  track$bacterial[is.na(track$bacterial)] <- 0
  write_table(track, file.path(stage, "read_tracking.tsv"))
  write_table(pct, file.path(stage, "phylum_percent.tsv"))
  write_table(runs, file.path(stage, "used_runs.tsv"))
  logs <- list.files(work, pattern = "_cutadapt\\.log$", full.names = TRUE)
  if (length(logs)) {
    dir.create(file.path(stage, "trimming_logs"))
    stopifnot(all(file.copy(logs, file.path(stage, "trimming_logs"))))
  }
  saveRDS(counts, file.path(stage, "asv_counts.rds"))
  saveRDS(tax, file.path(stage, "taxonomy.rds"))
  writeLines(capture.output(sessionInfo()), file.path(stage, "session_info.txt"))
  saveRDS(signature, file.path(stage, "complete.rds"))
  if (!file.rename(stage, final)) stop("Could not publish completed output folder: ", stage)
  invisible(final)
}

combine_studies <- function(projects = read_table("studies.tsv")$project, input = "results/full",
                            output = "results/combined") {
  if (!length(projects) || anyDuplicated(projects) || any(!grepl("^PRJ[A-Z]+[0-9]+$", projects))) stop("Choose valid unique project IDs")
  files <- file.path(input, projects, "phylum_percent.tsv")
  complete <- file.exists(file.path(input, projects, "complete.rds")) & file.exists(files)
  if (!all(complete)) stop("Missing completed results for ", paste(projects[!complete], collapse = ", "))
  records <- lapply(file.path(input, projects, "complete.rds"), readRDS)
  if (!all(vapply(records, function(x) identical(x$ref, reference_md5) && identical(x$scope, "full"), logical(1)))) {
    stop("Use completed full study results from this workflow with the same taxonomy reference")
  }
  data <- do.call(rbind, lapply(seq_along(files), function(i) {
    d <- read_table(files[i])
    if (!nrow(d) || any(d$project != projects[i])) stop("Unexpected project in ", files[i])
    d$percent <- as.numeric(d$percent)
    stopifnot(all(is.finite(d$percent) & d$percent >= 0), !anyDuplicated(d[c("sample", "phylum")]),
              all(abs(tapply(d$percent, d$sample, sum) - 100) < 1e-5))
    d
  }))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  write_table(data, file.path(output, "sample_phylum_percent.tsv"))
  n <- aggregate(sample ~ project, unique(data[c("project", "sample")]), length)
  mean_data <- aggregate(percent ~ project + phylum, data, sum)
  mean_data <- merge(expand.grid(project = projects, phylum = unique(data$phylum), stringsAsFactors = FALSE),
                     mean_data, all.x = TRUE)
  mean_data$percent[is.na(mean_data$percent)] <- 0
  mean_data$n_samples <- n$sample[match(mean_data$project, n$project)]
  mean_data$percent <- mean_data$percent / mean_data$n_samples
  write_table(mean_data, file.path(output, "study_mean_phylum_percent.tsv"))
  top <- head(names(sort(tapply(mean_data$percent, mean_data$phylum, mean), decreasing = TRUE)), 10)
  mean_data$phylum[!mean_data$phylum %in% top] <- "Other phyla"
  plotdata <- aggregate(percent ~ project + phylum, mean_data, sum)
  plotdata$project <- paste0(plotdata$project, "\n(n = ", n$sample[match(plotdata$project, n$project)], ")")
  plot <- ggplot2::ggplot(plotdata, ggplot2::aes(project, percent, fill = phylum)) +
    ggplot2::geom_col() + ggplot2::labs(x = NULL, y = "Mean within-sample bacterial read percentage",
      fill = "Phylum", title = "Gut bacterial composition by study",
      subtitle = "Descriptive comparison across studies and sequencing protocols") + ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(output, "study_phylum_composition.png"), plot, width = 11, height = 6, dpi = 200)
  invisible(mean_data)
}

run_team <- function(owner = "all", mode = c("inventory", "preview", "trial", "pilot", "full"),
                     max_download_gb = 2) {
  mode <- match.arg(mode)
  if (mode == "inventory") return(inventory())
  data <- load_inputs()
  keep <- owner == "all" | tolower(data$studies$owner) == tolower(owner) | data$studies$project == owner
  studies <- data$studies[keep, , drop = FALSE]
  if (!nrow(studies)) stop("Choose all, Vansh, Bhavya, Krupa, or a listed project ID")
  if (mode == "full" && any(studies$status != "ready")) {
    stop("Full run pending study review for ", paste(studies$project[studies$status != "ready"], collapse = ", "),
         ". See START_HERE.md and study_audit.md.")
  }
  selected <- lapply(seq_len(nrow(studies)), function(i) {
    runs <- runs_for_study(studies[i, , drop = FALSE], data$runs)
    if (mode %in% c("preview", "trial") && nrow(runs) > 3L) {
      runs <- runs[unique(round(seq(1, nrow(runs), length.out = 3))), ]
    }
    runs
  })
  todo <- seq_len(nrow(studies))
  if (mode == "full") {
    complete <- vapply(todo, function(i) full_result_available(studies[i, , drop = FALSE], selected[[i]]), logical(1))
    for (project in studies$project[complete]) message("Using completed results for ", project)
    todo <- todo[!complete]
    if (!length(todo)) {
      if (owner == "all") combine_studies()
      return(invisible(NULL))
    }
    existing <- studies$project[todo][dir.exists(file.path("results", "full", studies$project[todo]))]
    if (length(existing)) stop("Check existing full result files or settings for ", paste(existing, collapse = ", "))
  }
  if (!requireNamespace("dada2", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install DADA2 and ggplot2 as described in START_HERE.md")
  }
  if (packageVersion("dada2") != "1.38.0") stop("Use DADA2 1.38.0 for this workflow")
  if (mode == "pilot") {
    runs <- data$runs[match(pilot_ids, data$runs$run), ]
    stopifnot(identical(runs$run, pilot_ids))
    cfg <- data$studies[data$studies$project == "PRJNA687346", ]
    local <- file.path("pilot_input", paste0(pilot_ids, ".trimmed.fastq.gz"))
    stopifnot(all(file.exists(local)))
    fetch_file(reference_url, reference_file, reference_md5)
    return(study_results(cfg, runs, "results/pilot", local_trimmed = local))
  }
  gb <- sum(vapply(selected[todo], function(x) sum(as.numeric(c(x$bytes_f, x$bytes_r)), na.rm = TRUE), numeric(1))) / 1e9
  if (!is.finite(max_download_gb) || max_download_gb <= 0 || gb > max_download_gb) {
    stop(sprintf("Selected FASTQs total %.2f GB. Set max_download_gb above that if you have disk space for reads and intermediate files.", gb))
  }
  if (mode %in% c("trial", "full")) fetch_file(reference_url, reference_file, reference_md5)
  for (i in todo) {
    if (!nrow(selected[[i]])) {
      if (mode == "full") stop("No eligible runs for ", studies$project[i])
      message("Skipping ", studies$project[i], " until run selection is confirmed.")
      next
    }
    if (mode == "preview") {
      folder <- file.path("results", "preview", studies$project[i])
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      preview_project(selected[[i]], folder)
    } else study_results(studies[i, ], selected[[i]], file.path("results", mode))
  }
  if (mode == "full" && owner == "all") combine_studies()
  invisible(NULL)
}
