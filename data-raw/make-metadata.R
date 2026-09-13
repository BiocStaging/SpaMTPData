#!/usr/bin/env Rscript

# Delegate to the shipped recipe so the two entry points cannot drift.
argument <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script <- normalizePath(sub("^--file=", "", argument), mustWork = TRUE)
recipe <- normalizePath(file.path(dirname(script), "..", "inst", "scripts",
    "make-metadata.R"), mustWork = TRUE)
status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(recipe), shQuote(commandArgs(trailingOnly = TRUE))))
quit(save = "no", status = status)
