#!/usr/bin/env Rscript

# Generate candidate registry/Hub tables only after publication and fresh-cache reads.
# Usage: register-native-release.R SHARED_STAGING PUBLISHED_RECORD_ID OUTPUT
# Does not modify installed/package registries or perform Hub ingestion.
runNativeRegistration <- function(args) {
    if (length(args) != 3L) stop(paste("Usage: register-native-release.R",
        "SHARED_STAGING PUBLISHED_RECORD_ID OUTPUT"), call. = FALSE)
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) stop("Install SpatialExperiment.")
    if (!grepl("^[0-9]+$", args[[2L]])) stop("Use a numeric published record ID.")
    script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(),
        value = TRUE)[1L]), mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
    source(file.path(dirname(script), "native-registration-utils.R"), local = TRUE)
    output <- args[[3L]]
    nativeRequire(!dir.exists(output) || !length(list.files(output, all.files = TRUE,
        no.. = TRUE)), "Output directory is not empty; choose a new directory.")
    config <- jsonlite::read_json(file.path(args[[1L]], "release-config.json"))
    native <- utils::read.csv(file.path(args[[1L]], "native_resource_manifest.csv"),
        stringsAsFactors = FALSE)
    registry <- utils::read.csv(file.path(root, "inst", "manifest", "resource_manifest.csv"),
        stringsAsFactors = FALSE)
    metadata <- utils::read.csv(file.path(root, "inst", "extdata", "metadata.csv"),
        stringsAsFactors = FALSE)
    cache <- tempfile("spamtpdata-publication-check-")
    dir.create(cache)
    on.exit(unlink(cache, recursive = TRUE), add = TRUE)
    oldTimeout <- options(timeout = max(1800, getOption("timeout", 60)))
    on.exit(options(oldTimeout), add = TRUE)
    json <- file.path(cache, "record.json")
    utils::download.file(paste0("https://zenodo.org/api/records/", args[[2L]]),
        json, mode = "wb", quiet = TRUE)
    record <- jsonlite::read_json(json)
    nativeRequire(identical(as.character(record$id), args[[2L]]), "Record ID mismatch.")
    nativeRequire(identical(record$metadata$version, config$collection_version) &&
        identical(as.character(unique(native$version)), config$native_version),
        "Published collection or native resource version differs from the staged plan.")
    candidate <- buildNativeRegistration(registry, metadata, native, record,
        config$concept_record_id, config$base_record_id)
    for (i in seq_len(nrow(native))) {
        path <- file.path(cache, native$file_name[i])
        utils::download.file(paste0("https://zenodo.org/api/records/", args[[2L]],
            "/files/", native$file_name[i], "/content"), path, mode = "wb", quiet = TRUE)
        nativeRequire(file.info(path)$size == native$bytes[i] &&
            unname(tools::md5sum(path)) == native$md5[i], "Fresh download failed checksum.")
        object <- readRDS(path)
        nativeRequire(inherits(object, "SpatialExperiment") &&
            isTRUE(methods::validObject(object, test = TRUE)),
            "Fresh download is not a valid SpatialExperiment.")
        nativeRequire(!any(c("Seurat", "SeuratObject") %in% loadedNamespaces()),
            "Native reading unexpectedly loaded Seurat; aborting registration.")
        rm(object)
        invisible(gc())
    }
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(candidate$registry, file.path(output, "resource_manifest.csv"),
        row.names = FALSE, na = "")
    utils::write.csv(candidate$metadata, file.path(output, "metadata.csv"),
        row.names = FALSE, na = "")
    jsonlite::write_json(list(record_id = args[[2L]],
        resource_version = native$version[1L], verified_native_files = nrow(native),
        complete_snapshot_resources = sum(candidate$registry$version == native$version[1L]),
        hub_ingestion = "not performed", package_registries = "unchanged"),
        file.path(output, "registration-validation.json"), pretty = TRUE, auto_unbox = TRUE)
    message("Fresh downloads verified. Candidate tables written to ", normalizePath(output),
        "; review before replacing package tables and requesting Hub ingestion.")
}
runNativeRegistration(commandArgs(trailingOnly = TRUE))
