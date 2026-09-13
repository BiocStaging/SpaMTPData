#!/usr/bin/env Rscript

if (!requireNamespace("ExperimentHubData", quietly = TRUE)) {
    stop("Install ExperimentHubData before validating Hub metadata.")
}
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1L]
script_path <- sub("^--file=", "", script_arg)
package_root <- normalizePath(
    file.path(dirname(script_path), "..", ".."), mustWork = TRUE
)
manifest <- utils::read.csv(file.path(package_root, "inst", "manifest",
    "resource_manifest.csv"), stringsAsFactors = FALSE)
metadata <- utils::read.csv(file.path(package_root, "inst", "extdata",
    "metadata.csv"), stringsAsFactors = FALSE)
index <- match(manifest$title, metadata$Title)
stopifnot(!anyNA(index), !anyDuplicated(metadata$Title),
    nrow(manifest) == nrow(metadata),
    identical(manifest$r_data_class, metadata$RDataClass[index]),
    identical(manifest$dispatch_class, metadata$DispatchClass[index]),
    identical(manifest$rdata_path, metadata$RDataPath[index]),
    identical(manifest$location_prefix, metadata$Location_Prefix[index]))
ExperimentHubData::makeExperimentHubMetadata(package_root, fileName = "metadata.csv")
