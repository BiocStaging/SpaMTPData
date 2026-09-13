#!/usr/bin/env Rscript

# One-time preparation; NOT a runtime dependency of the data package.
# Usage: Rscript inst/scripts/prepare-native-resources.R SOURCES OUTPUT [1.1.0]
# SOURCES: directory of official archives, or CSV with resource/version/path.
# OUTPUT: staging root; native files are written under OUTPUT/VERSION.
# Requires an installed SpaMTP >= 0.99.4 and optional Seurat for old image
# class definitions. Nothing is uploaded and published registries are untouched.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
    stop("Usage: prepare-native-resources.R SOURCES OUTPUT [NATIVE_VERSION]")
}
version <- if (length(args) == 3L) args[[3L]] else "1.1.0"
if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)) stop("Use a numeric x.y.z version.")
for (package in c("SpaMTP", "Seurat", "SeuratObject")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        stop("One-time resource preparation requires ", package,
             "; native output files do not depend on Seurat.")
    }
}
if (utils::packageVersion("SpaMTP") < "0.99.4") stop("Use SpaMTP >= 0.99.4.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
script <- normalizePath(script, mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), "..", ".."))
helper <- file.path(root, "inst", "scripts", "native-resource-utils.R")
planFile <- file.path(root, "inst", "manifest", "native_conversion_plan.csv")
source(helper, local = TRUE)
plan <- utils::read.csv(planFile, stringsAsFactors = FALSE)
registry <- utils::read.csv(file.path(root, "inst", "manifest", "resource_manifest.csv"),
    stringsAsFactors = FALSE)
registry <- registry[registry$r_data_class == "Seurat", , drop = FALSE]
stopifnot(nrow(plan) == 7L, !anyDuplicated(plan$resource),
    setequal(plan$resource, registry$resource))
registry <- registry[match(plan$resource, registry$resource), , drop = FALSE]
if (dir.exists(args[[1L]])) {
    sourceRoot <- normalizePath(args[[1L]])
    paths <- vapply(seq_len(nrow(registry)), function(i) {
        candidates <- c(file.path(sourceRoot, registry$version[i], registry$file_name[i]),
            file.path(sourceRoot, registry$file_name[i]))
        existing <- candidates[file.exists(candidates)]
        if (!length(existing)) stop("Missing source archive: ", registry$file_name[i])
        existing[[1L]]
    }, character(1))
} else {
    input <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)
    stopifnot(all(c("resource", "version", "path") %in% names(input)),
        !anyDuplicated(paste(input$resource, input$version)))
    index <- match(paste(registry$resource, registry$version), paste(input$resource, input$version))
    stopifnot(!anyNA(index))
    paths <- input$path[index]
}
stopifnot(all(file.exists(paths)),
    identical(unname(tools::md5sum(paths)), registry$md5),
    all(file.info(paths)$size == registry$bytes))
destination <- file.path(args[[2L]], version)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
destination <- normalizePath(destination, mustWork = TRUE)
recipe <- unname(tools::md5sum(c(script, helper, planFile)))
names(recipe) <- c("script", "helper", "plan")
versions <- vapply(c("SpaMTP", "Seurat", "SeuratObject", "SpatialExperiment",
    "SingleCellExperiment", "SummarizedExperiment", "Matrix"), function(package) {
    as.character(utils::packageVersion(package))
}, character(1))
records <- vector("list", nrow(plan))
for (i in seq_len(nrow(plan))) {
    specification <- plan[i, , drop = FALSE]
    row <- registry[i, , drop = FALSE]
    cat("RESOURCE", row$resource, "\n")
    original <- readRDS(paths[i])
    filename <- paste0(row$resource, "_spe.rds")
    target <- file.path(destination, filename)
    if (file.exists(target)) {
        object <- readRDS(target)
        provenance <- S4Vectors::metadata(object)$native_resource
        if (!identical(provenance$source_md5, row$md5) ||
            !identical(provenance$recipe_md5, recipe) ||
            !identical(provenance$package_versions, versions) ||
            !identical(provenance$native_version, version)) {
            stop("An existing native file has different provenance. Use a new staging directory: ", target)
        }
    } else {
        object <- prepareNativeResource(original, specification)
        S4Vectors::metadata(object)$native_resource <- list(
            resource = row$resource, native_version = version,
            source_version = row$version, source_file = row$file_name,
            source_class = row$r_data_class, source_url = row$source_url,
            source_md5 = row$md5, source_bytes = row$bytes,
            recipe_md5 = recipe, package_versions = versions,
            r_version = as.character(getRversion()), status = "locally prepared; not uploaded")
        validateNativeResource(object, original, specification)
        temporary <- tempfile(paste0(row$resource, "-"), tmpdir = destination, fileext = ".part")
        saveRDS(object, temporary, compress = "gzip", version = 3L)
        restored <- readRDS(temporary)
        stopifnot(identical(object, restored))
        validateNativeResource(restored, original, specification)
        if (!file.rename(temporary, target)) stop("Cannot move verified native file into place: ", target)
    }
    validateNativeResource(object, original, specification)
    records[[i]] <- data.frame(resource = row$resource, version = version,
        file_name = filename, r_data_class = class(object)[1L],
        rows = nrow(object), columns = ncol(object),
        primary_assay = specification$primary_assay,
        alternative_assays = paste(SingleCellExperiment::altExpNames(object), collapse = ";"),
        expression_layers = paste(SummarizedExperiment::assayNames(object), collapse = ";"),
        source_image = specification$image,
        source_version = row$version, source_file = row$file_name,
        source_md5 = row$md5, source_bytes = row$bytes, source_url = row$source_url,
        bytes = file.info(target)$size, md5 = unname(tools::md5sum(target)),
        publication_status = "local", stringsAsFactors = FALSE)
    cat("SAVED", filename, nrow(object), "x", ncol(object),
        "ALT", paste(SingleCellExperiment::altExpNames(object), collapse = ","),
        "BYTES", file.info(target)$size, "\n")
    rm(original, object)
    invisible(gc())
}
manifest <- do.call(rbind, records)
utils::write.csv(manifest, file.path(destination, "native_resource_manifest.csv"), row.names = FALSE)
cat("Prepared and verified", nrow(manifest), "native RDS files in", destination, "\n")
