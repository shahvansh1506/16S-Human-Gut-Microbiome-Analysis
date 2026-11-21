# Open gut16s-team.Rproj. Put the six completed study folders under
# results/full. Run source("make_final_report.R") to generate the report.
source("workflow.R")

make_final_report <- function(input = "results/full", output = "results/combined") {
  studies <- read_table("studies.tsv")
  if (nrow(studies) != 6L || anyDuplicated(studies$project)) {
    stop("This report expects six distinct studies in studies.tsv.")
  }
  files <- c("complete.rds", "read_tracking.tsv", "used_runs.tsv", "phylum_percent.tsv")
  paths <- file.path(input, studies$project)
  ready <- vapply(paths, function(path) all(file.exists(file.path(path, files))), logical(1))
  if (!all(ready)) stop("Missing study results: ", paste(studies$project[!ready], collapse = ", "))

  # Verify each folder's recorded study and sample selection before combining.
  records <- lapply(paths, function(path) readRDS(file.path(path, "complete.rds")))
  audit <- lapply(seq_len(nrow(studies)), function(i) {
    project <- studies$project[i]
    folder <- paths[i]
    record <- records[[i]]
    if (!identical(record$scope, "full") || !identical(record$study$project, project)) {
      stop("Expected full output for ", project)
    }
    used <- read_table(file.path(folder, "used_runs.tsv"))
    track <- read_table(file.path(folder, "read_tracking.tsv"))
    comp <- read_table(file.path(folder, "phylum_percent.tsv"))
    if (!setequal(record$runs$run, used$run) || !setequal(used$run, track$run) ||
        anyDuplicated(track$run) || !all(comp$run %in% used$run)) {
      stop("Run tracking and recorded inputs disagree for ", project)
    }
    incoming <- as.numeric(track$input)
    retained <- as.numeric(track$nonchim)
    if (any(!is.finite(incoming) | !is.finite(retained) | incoming < 0 | retained < 0 | retained > incoming)) {
      stop("Invalid read tracking counts for ", project)
    }
    data.frame(study = project, owner = studies$owner[i], analyzed_runs = nrow(used),
               samples_in_figure = length(unique(comp$sample)),
               median_nonchim_percent = round(median(100 * retained[incoming > 0] / incoming[incoming > 0]), 1),
               region = record$study$region, layout = paste(unique(used$layout), collapse = "/"))
  })
  summary <- do.call(rbind, audit)
  combine_studies(studies$project, input, output)
  write_table(summary, file.path(output, "study_read_summary.tsv"))
  means <- read_table(file.path(output, "study_mean_phylum_percent.tsv"))
  means$percent <- as.numeric(means$percent)
  top <- head(names(sort(tapply(means$percent, means$phylum, mean), decreasing = TRUE)), 10L)
  means$phylum[!means$phylum %in% top] <- "Other phyla"
  plot_data <- aggregate(percent ~ project + phylum, means, sum)
  label <- setNames(paste0(summary$study, "\n(n = ", summary$samples_in_figure, ")"), summary$study)
  figure <- ggplot2::ggplot(plot_data, ggplot2::aes(project, percent, fill = phylum)) +
    ggplot2::geom_col() + ggplot2::scale_x_discrete(labels = label) +
    ggplot2::labs(title = "Bacterial composition across six gut microbiome studies",
                  x = NULL, y = "Mean bacterial read percentage", fill = "Phylum") +
    ggplot2::theme_minimal(base_size = 10)

  grDevices::pdf(file.path(output, "six_study_results_report.pdf"), width = 11, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(2, 1, 2, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  line_at <- function(text, y, x = 0.04, size = 1, bold = FALSE) {
    graphics::text(x, y, text, adj = c(0, 1), cex = size, font = if (bold) 2 else 1, col = "#17354A")
  }
  line_at("Gut microbiome | Six-study results", 0.94, size = 2, bold = TRUE)
  line_at("Vansh Shah, Bhavya and Krupa", 0.86, size = 1.1)
  line_at("Workflow and retained samples", 0.77, bold = TRUE)
  line_at("One R workflow applied study-specific trimming and DADA2 settings and assigned taxonomy with SILVA 138.2.", 0.73, size = 0.92)
  line_at("Bacterial phylum percentages were calculated within each sample and averaged with equal sample weights.", 0.695, size = 0.92)
  columns <- c(0.04, 0.24, 0.39, 0.59, 0.84)
  headings <- c("Study", "Region", "Runs processed", "Samples in figure", "Median retained*")
  for (j in seq_along(columns)) line_at(headings[j], 0.615, columns[j], 0.9, TRUE)
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    values <- c(row$study, row$region, row$analyzed_runs, row$samples_in_figure, paste0(row$median_nonchim_percent, "%"))
    for (j in seq_along(columns)) line_at(values[j], 0.575-i*0.043, columns[j], 0.9)
  }
  line_at("* Median nonchimera reads divided by reads entering DADA2, after primer trimming when used.", 0.29, size = 0.83)
  line_at("Interpretation", 0.23, bold = TRUE)
  line_at("This chart compares observed bacterial composition. Cohort and laboratory differences can influence the percentages.", 0.19, size = 0.90)
  line_at("Review the recorded inclusion decisions, archive notes and read-tracking files alongside these results.", 0.15, size = 0.90)
  noqual <- summary$study[vapply(records, function(record) identical(record$study$error_model, "constant_quality"), logical(1))]
  if (length(noqual)) {
    line_at(paste0("Models without quality-score dependence: ", paste(noqual, collapse = ", "), "."), 0.11, size = 0.90)
  }
  print(figure)
  message("Results report saved in ", output)
  invisible(summary)
}

make_final_report()
