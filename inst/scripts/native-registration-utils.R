# Maintainer-only registration helpers; no dependency on SpaMTPdb or Seurat.

nativeRequire <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
}

buildNativeRegistration <- function(registry, metadata, native, record,
                                    conceptId, baseRecordId) {
    nativeRequire(isTRUE(record$submitted) && identical(record$status, "published") &&
        identical(record$metadata$access_right, "open"),
        "Only a published, publicly accessible record can be registered.")
    id <- as.character(record$id)
    nativeRequire(length(id) == 1L && grepl("^[0-9]+$", id) &&
        id != as.character(baseRecordId) &&
        identical(as.character(record$conceptrecid), as.character(conceptId)),
        "Record is not a new version in the shared Zenodo family.")
    required <- c("resource", "version", "file_name", "r_data_class", "bytes", "md5",
        "source_version", "source_file", "source_md5", "source_bytes", "source_url")
    nativeRequire(all(required %in% names(native)) && nrow(native) > 0L &&
        !anyNA(native[, required]) && !anyDuplicated(native$resource) &&
        !anyDuplicated(native$file_name), "Incomplete or duplicated native manifest.")
    nativeRequire(length(unique(native$version)) == 1L &&
        length(unique(native$source_version)) == 1L &&
        utils::compareVersion(native$version[1L], native$source_version[1L]) > 0L,
        "Select one new native version and one source registry version.")
    source <- registry[registry$version == native$source_version[1L], , drop = FALSE]
    nativeRequire(nrow(source) > 0L && !anyDuplicated(source$resource) &&
        setequal(native$resource, source$resource[source$r_data_class == "Seurat"]),
        "Native replacement set does not cover every historical Seurat resource.")
    index <- match(native$resource, source$resource)
    nativeRequire(identical(native$source_md5, source$md5[index]) &&
        identical(native$source_file, source$file_name[index]) &&
        identical(native$source_url, source$source_url[index]) &&
        all(native$source_bytes == source$bytes[index]) &&
        all(native$r_data_class == "SpatialExperiment") &&
        all(grepl("^[A-Za-z0-9][A-Za-z0-9._-]*[.]rds$", native$file_name)) &&
        all(grepl("^[0-9a-f]{32}$", native$md5)) && all(native$bytes > 0),
        "Native provenance, filename, checksum or container class is invalid.")
    remote <- do.call(rbind, lapply(record$files, function(x) {
        data.frame(file_name = x$key, bytes = x$size,
            md5 = sub("^md5:", "", x$checksum), stringsAsFactors = FALSE)
    }))
    nativeRequire(!is.null(remote) && !anyDuplicated(remote$file_name),
        "Published file list is missing or duplicated.")
    fileIndex <- match(native$file_name, remote$file_name)
    nativeRequire(!anyNA(fileIndex) && identical(native$md5, remote$md5[fileIndex]) &&
        all(native$bytes == remote$bytes[fileIndex]),
        "Published native file sizes or MD5 do not match the prepared files.")
    metaIndex <- match(source$title, metadata$Title)
    nativeRequire(!anyNA(metaIndex) && !anyDuplicated(metadata$Title),
        "Source Hub metadata is missing or duplicated.")
    # A release is a complete snapshot: carry forward ALL unchanged resources.
    nextRegistry <- source
    nextRegistry$version <- native$version[1L]
    nextRegistry$title <- paste0("SpaMTPData_", source$resource, "_", native$version[1L])
    for (column in c("file_name", "r_data_class", "bytes", "md5")) {
        nextRegistry[index, column] <- native[[column]]
    }
    prefix <- paste0("https://zenodo.org/api/records/", id, "/files/")
    nextRegistry$location_prefix[index] <- prefix
    nextRegistry$rdata_path[index] <- paste0(native$file_name, "/content")
    nextRegistry$source_url[index] <- paste0("https://zenodo.org/records/", id)
    nextRegistry$dispatch_class[index] <- "Rds"
    nextMetadata <- metadata[metaIndex, , drop = FALSE]
    rownames(nextMetadata) <- NULL
    nextMetadata$Title <- nextRegistry$title
    nextMetadata$RDataClass <- nextRegistry$r_data_class
    nextMetadata$DispatchClass <- nextRegistry$dispatch_class
    nextMetadata$Location_Prefix <- nextRegistry$location_prefix
    nextMetadata$RDataPath <- nextRegistry$rdata_path
    nextMetadata$Description[index] <- paste(nextMetadata$Description[index],
        "Native SpatialExperiment retaining selected assays, annotations and centroids;",
        "optical, polygon and Seurat graph/reduction payloads are omitted.")
    nextMetadata$Tags[index] <- paste0(nextMetadata$Tags[index], ":SpatialExperiment:Native")
    # SourceUrl/SourceVersion retain upstream provenance, not the download location.
    nextMetadata$SourceUrl[index] <- native$source_url
    nextMetadata$SourceVersion[index] <- native$source_version
    existing <- registry[registry$version == native$version[1L], , drop = FALSE]
    if (nrow(existing)) {
        existing <- existing[match(nextRegistry$resource, existing$resource), , drop = FALSE]
        rownames(existing) <- rownames(nextRegistry) <- NULL
        existingMetadata <- metadata[match(nextMetadata$Title, metadata$Title), , drop = FALSE]
        rownames(existingMetadata) <- NULL
        nativeRequire(identical(existing, nextRegistry) &&
            identical(existingMetadata, nextMetadata),
            "That resource version already exists with different contents or metadata.")
        return(list(registry = registry, metadata = metadata))
    }
    nativeRequire(!any(nextMetadata$Title %in% metadata$Title),
        "A Hub title for this release already exists.")
    list(registry = rbind(registry, nextRegistry), metadata = rbind(metadata, nextMetadata))
}
