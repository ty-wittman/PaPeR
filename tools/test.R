#!/usr/bin/env Rscript
# Run the PaPeR test suite against an optimized install of the working tree.
#
#   Rscript tools/test.R                     # quick tier (about a minute)
#   Rscript tools/test.R --tier=standard     # all supported tests
#   Rscript tools/test.R --tier=extended     # everything
#   Rscript tools/test.R --files=test-public-api.R,test-stats.R
#   Rscript tools/test.R --lib=/path/to/lib  # test an existing install
#
# Options: --jobs=N (parallel test files, default 4); --out=FILE writes one
# row per file (seconds, passed, failed, errors, skipped). Each test file runs in
# the PaPeR namespace with BLAS/OpenMP threads at one. Exit status is nonzero
# when any expectation fails or errors.

args <- commandArgs(TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
repo <- dirname(dirname(script))
pkg <- repo
test_dir <- file.path(pkg, "tests", "testthat")

# All supported test files currently fit in the quick tier. Keep the tier
# interface so slower supported files can be assigned separately when needed.
standard_only <- character()
extended_only <- character()
all_files <- sort(list.files(test_dir, "^test-.*[.]R$"))
tier <- opt("tier", "quick")
if (!tier %in% c("quick", "standard", "extended")) stop("--tier must be quick, standard or extended")
files <- switch(tier,
  quick = setdiff(all_files, c(standard_only, extended_only)),
  standard = setdiff(all_files, extended_only),
  extended = all_files)
if (!is.null(opt("files"))) {
  files <- strsplit(opt("files"), ",", fixed = TRUE)[[1]]; tier <- "selected"
}
missing <- setdiff(files, all_files)
if (length(missing)) stop("unknown test file(s): ", paste(missing, collapse = ", "))
jobs <- as.integer(opt("jobs", "4"))

lib <- opt("lib")
if (is.null(lib)) {
  lib <- file.path(tempdir(), "pap-test-lib"); dir.create(lib, showWarnings = FALSE)
  cat("installing", pkg, "\n")
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-test-load", paste0("--library=", shQuote(lib)), shQuote(pkg)),
    stdout = file.path(lib, "install.log"), stderr = file.path(lib, "install.log"),
    env = "MAKEFLAGS=-j4")
  if (status != 0) stop("install failed; see ", file.path(lib, "install.log"))
}
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", TESTTHAT_CPUS = "1")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(testthat); library(PaPeR, lib.loc = lib) })
setwd(repo)

run_one <- function(f) {
  started <- proc.time()[["elapsed"]]
  r <- tryCatch(as.data.frame(test_file(file.path(test_dir, f), reporter = "silent",
                                        stop_on_failure = FALSE, package = "PaPeR",
                                        load_package = "none")),
                error = function(e) data.frame(test = conditionMessage(e), failed = 0L,
                                               error = TRUE, passed = 0L, skipped = 0L))
  data.frame(file = f, seconds = round(proc.time()[["elapsed"]] - started, 1),
             passed = sum(r$passed), failed = sum(r$failed), errors = sum(r$error),
             skipped = sum(r$skipped),
             problems = paste(r$test[r$failed > 0 | r$error], collapse = "; "))
}
cat(sprintf("running %d test files (%s, %d jobs)\n", length(files), tier, jobs))
started <- proc.time()[["elapsed"]]
rows <- do.call(rbind, parallel::mclapply(files, run_one, mc.cores = jobs, mc.preschedule = FALSE))
if (!is.null(opt("out"))) write.table(rows, opt("out"), sep = "\t", quote = FALSE, row.names = FALSE)
slow <- head(rows[order(-rows$seconds), ], 3)
cat("slowest:", paste0(slow$file, " (", slow$seconds, " s)", collapse = ", "), "\n")
bad <- rows[rows$failed > 0 | rows$errors > 0, , drop = FALSE]
cat(sprintf("%d files, %d passed, %d failed, %d errors, %d skipped in %.0f s\n",
            nrow(rows), sum(rows$passed), sum(rows$failed), sum(rows$errors),
            sum(rows$skipped), proc.time()[["elapsed"]] - started))
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) cat("FAIL", bad$file[i], ":", bad$problems[i], "\n")
  quit(status = 1)
}
