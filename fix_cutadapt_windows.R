# Load this after workflow.R in the same Windows RStudio session.
cutadapt_runner <- local({
  original_runner <- cutadapt_runner

  function() {
    wsl <- unname(Sys.which("wsl.exe"))
    if (!nzchar(wsl)) wsl <- file.path(Sys.getenv("SystemRoot"), "System32", "wsl.exe")
    if (.Platform$OS.type != "windows" || !file.exists(wsl)) return(original_runner())

    # WSL sometimes writes its distribution list as UTF-16LE on Windows.
    listing <- tempfile("gut16s-wsl-")
    on.exit(unlink(listing), add = TRUE)
    status <- system2(wsl, c("--list", "--quiet"), stdout = listing, stderr = FALSE)
    if (status != 0L || !file.exists(listing)) {
      stop("Windows could not list Ubuntu installations. Run wsl -l -v in PowerShell.")
    }
    bytes <- readBin(listing, "raw", n = file.size(listing))
    encoding <- if (any(bytes == as.raw(0))) "UTF-16LE" else "UTF-8"
    names <- iconv(list(bytes), from = encoding, to = "UTF-8")
    names <- trimws(strsplit(sub("^\ufeff", "", names), "\n", fixed = TRUE)[[1L]])
    distros <- grep("^Ubuntu", names, ignore.case = TRUE, value = TRUE)
    if (!length(distros)) {
      stop("Windows did not list any Ubuntu installations. Run wsl -l -v in PowerShell.")
    }

    project <- getwd()
    invoke_wsl <- function(args) {
      needs_quotes <- grepl("[[:space:]\"]", args)
      args[needs_quotes] <- shQuote(args[needs_quotes])
      suppressWarnings(system2(wsl, args, stdout = TRUE, stderr = TRUE))
    }
    details <- character()
    for (distro in distros) {
      prefix <- c("-d", distro, "--cd", project, "--exec")
      home <- invoke_wsl(c(prefix, "/usr/bin/printenv", "HOME"))
      if (!is.null(attr(home, "status")) || length(home) != 1L || !startsWith(home, "/")) {
        details <- c(details, paste(distro, paste(home, collapse = " ")))
        next
      }
      for (executable in c(file.path(home, "gut16s-cutadapt/bin/cutadapt"), "cutadapt")) {
        runner <- local({
          command <- c(prefix, executable)
          function(args) {
            args <- chartr("\\", "/", args)
            invoke_wsl(c(command, args))
          }
        })
        version <- runner("--version")
        if (is.null(attr(version, "status")) && any(trimws(version) == "4.7")) {
          message("RStudio found Cutadapt 4.7 in ", distro)
          return(runner)
        }
        details <- c(details, paste(distro, executable, paste(version, collapse = " ")))
      }
    }
    stop("RStudio could not start Cutadapt 4.7.\n", paste(details, collapse = "\n"))
  }
})
