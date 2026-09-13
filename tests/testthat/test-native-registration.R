nativeRegistrationFixture <- function(path) {
    dir.create(path)
    environment <- new.env(parent = baseenv())
    sys.source(system.file("scripts", "native-registration-utils.R", package = "SpaMTPData"),
        envir = environment)
    registry <- utils::read.csv(system.file("manifest", "resource_manifest.csv",
        package = "SpaMTPData"), stringsAsFactors = FALSE)
    metadata <- utils::read.csv(system.file("extdata", "metadata.csv",
        package = "SpaMTPData"), stringsAsFactors = FALSE)
    # Build the transition from its historical snapshot, not the live latest.
    registry <- registry[registry$version == "1.0.0", , drop = FALSE]
    metadata <- metadata[match(registry$title, metadata$Title), , drop = FALSE]
    rownames(registry) <- rownames(metadata) <- NULL
    old <- registry[registry$r_data_class == "Seurat", , drop = FALSE]
    native <- old[, c("resource", "version", "file_name", "r_data_class", "bytes", "md5")]
    native$version <- "1.1.0"
    native$file_name <- paste0(native$resource, "_spe.rds")
    native$r_data_class <- "SpatialExperiment"
    object <- spaMTPExampleData()
    for (name in native$file_name) saveRDS(object, file.path(path, name))
    native$bytes <- file.info(file.path(path, native$file_name))$size
    native$md5 <- unname(tools::md5sum(file.path(path, native$file_name)))
    native$source_version <- old$version
    native$source_file <- old$file_name
    native$source_md5 <- old$md5
    native$source_bytes <- old$bytes
    native$source_url <- old$source_url
    record <- list(id = 999, conceptrecid = "22045310", submitted = TRUE,
        status = "published", metadata = list(access_right = "open"),
        files = lapply(seq_len(nrow(native)), function(i) list(key = native$file_name[i],
            size = native$bytes[i], checksum = paste0("md5:", native$md5[i]))))
    list(h = environment, registry = registry, metadata = metadata, native = native,
        object = object, record = record)
}

test_that("native registration keeps the complete 18-resource snapshot and history", {
    path <- tempfile("native-registration-")
    f <- nativeRegistrationFixture(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    result <- f$h$buildNativeRegistration(f$registry, f$metadata, f$native,
        f$record, "22045310", "22045311")
    expect_equal(nrow(result$registry), 36L)
    expect_equal(nrow(result$metadata), 36L)
    # file.info() uses doubles while CSV sizes may be integers; preserve values.
    expect_equal(result$registry[seq_len(18L), ], f$registry)
    expect_identical(result$metadata[seq_len(18L), ], f$metadata)
    latest <- result$registry[result$registry$version == "1.1.0", ]
    expect_setequal(latest$resource, f$registry$resource)
    expect_equal(sum(latest$r_data_class == "SpatialExperiment"), 7L)
    expect_false(any(latest$r_data_class == "Seurat"))
    unchanged <- !f$registry$resource %in% f$native$resource
    for (column in c("location_prefix", "rdata_path", "md5", "bytes", "r_data_class")) {
        expect_equal(latest[[column]][unchanged], f$registry[[column]][unchanged])
    }
    native <- latest$resource %in% f$native$resource
    expect_true(all(latest$location_prefix[native] == "https://zenodo.org/api/records/999/files/"))
    again <- f$h$buildNativeRegistration(result$registry, result$metadata,
        f$native, f$record, "22045310", "22045311")
    expect_identical(again, result)
})

test_that("the existing API can read native resources without a database-package dependency", {
    path <- tempfile("native-registration-")
    f <- nativeRegistrationFixture(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    result <- f$h$buildNativeRegistration(f$registry, f$metadata, f$native,
        f$record, "22045310", "22045311")
    manifest <- result$registry
    manifest$requires_conversion <- manifest$r_data_class == "Seurat"
    local_mocked_bindings(.spamtpdata_manifest = function() manifest)
    expect_equal(nrow(spaMTPDataResources(version = "latest")), 18L)
    for (resource in f$native$resource) {
        observed <- spaMTPData(resource, local_dir = path, offline = TRUE)
        provenance <- S4Vectors::metadata(observed)$SpaMTPData
        expect_identical(provenance$resource, resource)
        expect_identical(provenance$version, "1.1.0")
        expect_identical(provenance$md5, f$native$md5[match(resource, f$native$resource)])
        child <- SingleCellExperiment::altExp(observed, "transcriptome")
        expect_identical(S4Vectors::metadata(child)$SpaMTPData, provenance)
        S4Vectors::metadata(child)$SpaMTPData <- NULL
        SingleCellExperiment::altExp(observed, "transcriptome") <- child
        S4Vectors::metadata(observed)$SpaMTPData <- NULL
        expect_identical(observed, f$object)
    }
    historical <- spaMTPData(f$native$resource[1L], version = "1.0.0", metadata = TRUE)
    expect_identical(historical$r_data_class, "Seurat")
    expect_identical(historical$md5, f$native$source_md5[1L])
    expect_false(grepl("SpaMTPdb|Seurat", utils::packageDescription("SpaMTPData", fields = "Imports")))
})

test_that("registration rejects drafts, other families, damaged files and partial migrations", {
    path <- tempfile("native-registration-")
    f <- nativeRegistrationFixture(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    build <- function(record = f$record, native = f$native) {
        f$h$buildNativeRegistration(f$registry, f$metadata, native, record,
            "22045310", "22045311")
    }
    record <- f$record
    record$submitted <- FALSE
    expect_error(build(record), "published")
    record <- f$record
    record$conceptrecid <- "123"
    expect_error(build(record), "family")
    record <- f$record
    record$id <- 22045311
    expect_error(build(record), "new version")
    record <- f$record
    record$files[[1L]]$checksum <- strrep("0", 32L)
    expect_error(build(record), "MD5")
    record <- f$record
    record$files <- record$files[-1L]
    expect_error(build(record), "MD5")
    expect_error(build(native = f$native[-1L, ]), "every historical")
    native <- f$native
    native$source_md5[1L] <- strrep("0", 32L)
    expect_error(build(native = native), "provenance")
    native <- f$native
    native$r_data_class[1L] <- "Seurat"
    expect_error(build(native = native), "container class")
})
