.spamtpdata_manifest <- function() {
    path <- system.file("manifest", "resource_manifest.csv", package = "SpaMTPData")
    if (!nzchar(path)) {
        stop("SpaMTPData resource manifest is unavailable.", call. = FALSE)
    }
    manifest <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    manifest$requires_conversion <- manifest$r_data_class == "Seurat"
    manifest
}

.spamtpdata_resolve_version <- function(manifest, version) {
    versions <- unique(as.character(manifest$version))
    if (is.null(version) || identical(version, "latest")) {
        return(tail(sort(package_version(versions)), 1L) |> as.character())
    }
    version <- as.character(version)[1L]
    if (!version %in% versions) {
        stop(
            "SpaMTPData version '", version, "' is unavailable. Available: ",
            paste(versions, collapse = ", "),
            call. = FALSE
        )
    }
    version
}

.spamtpdata_local_dir <- function(local_dir = NULL) {
    if (!is.null(local_dir)) {
        return(normalizePath(local_dir, mustWork = FALSE))
    }
    configured <- getOption("SpaMTPData.resource_dir", "")
    if (!nzchar(configured)) {
        configured <- Sys.getenv("SPAMTPDATA_RESOURCE_DIR", "")
    }
    if (!nzchar(configured)) NULL else normalizePath(configured, mustWork = FALSE)
}

.spamtpdata_local_file <- function(row, local_dir) {
    if (is.null(local_dir)) return(NULL)
    candidates <- unique(c(
        file.path(local_dir, row$version, row$file_name),
        file.path(local_dir, row$file_name),
        file.path(local_dir, row$resource, row$file_name),
        file.path(local_dir, paste0(row$resource, ".rds")),
        file.path(local_dir, paste0(row$resource, ".RDS"))
    ))
    found <- candidates[file.exists(candidates)]
    if (length(found)) found[[1L]] else NULL
}

.spamtpdata_read_local <- function(path, dispatch_class) {
    if (tolower(dispatch_class) %in% c("rds", "rda")) {
        if (tolower(dispatch_class) == "rds") return(readRDS(path))
        environment <- new.env(parent = emptyenv())
        loaded <- load(path, envir = environment)
        if (length(loaded) != 1L) {
            stop("Local Rda resource must contain exactly one object.", call. = FALSE)
        }
        return(environment[[loaded]])
    }
    normalizePath(path, mustWork = TRUE)
}

.spamtpdata_cache_dir <- function(cache_dir = NULL, create = TRUE) {
    if (is.null(cache_dir)) {
        cache_dir <- getOption("SpaMTPData.cache_dir", "")
    }
    if (!nzchar(cache_dir)) {
        cache_dir <- Sys.getenv("SPAMTPDATA_CACHE_DIR", "")
    }
    if (!nzchar(cache_dir)) {
        cache_dir <- tools::R_user_dir("SpaMTPData", which = "cache")
    }
    if (create) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    normalizePath(cache_dir, mustWork = create)
}

.spamtpdata_file_valid <- function(path, row) {
    if (!file.exists(path) || isTRUE(file.info(path)$isdir)) return(FALSE)
    expected_bytes <- suppressWarnings(as.numeric(row$bytes[[1L]]))
    if (length(expected_bytes) != 1L || !is.finite(expected_bytes) ||
        expected_bytes <= 0 || file.info(path)$size != expected_bytes) {
        return(FALSE)
    }
    expected_md5 <- tolower(as.character(row$md5[[1L]]))
    if (length(expected_md5) != 1L || is.na(expected_md5) ||
        !grepl("^[0-9a-f]{32}$", expected_md5)) return(FALSE)
    identical(tolower(unname(tools::md5sum(path))), expected_md5)
}

.spamtpdata_cached_file <- function(row, cache_dir = NULL) {
    root <- .spamtpdata_cache_dir(cache_dir, create = FALSE)
    paths <- c(file.path(root, row$version, row$file_name),
               file.path(root, row$file_name))
    valid <- vapply(paths, .spamtpdata_file_valid, logical(1), row = row)
    if (any(valid)) paths[which(valid)[1L]] else NULL
}

.spamtpdata_verify_local <- function(path, row, verify) {
    if (verify && !.spamtpdata_file_valid(path, row)) {
        stop("Local resource '", row$resource, "' failed its size or MD5 check. ",
             "Use verify = FALSE only for intentional development fixtures, ",
             "not to label modified data as an official release.", call. = FALSE)
    }
    path
}

.spamtpdata_check_class <- function(value, row) {
    if (identical(row$dispatch_class[[1L]], "FilePath")) {
        valid <- is.character(value) && length(value) == 1L &&
            !is.na(value) && file.exists(value)
    } else {
        expected <- strsplit(row$r_data_class[[1L]], "/", fixed = TRUE)[[1L]]
        valid <- any(vapply(expected, function(type) inherits(value, type), logical(1)))
    }
    if (!valid) {
        stop("Resource '", row$resource, "' has class ",
             paste(class(value), collapse = "/"), "; expected ", row$r_data_class,
             ".", call. = FALSE)
    }
    value
}

.spamtpdata_download <- function(row, cache_dir = NULL, timeout = 1800,
                                 retries = 3L) {
    cache_dir <- .spamtpdata_cache_dir(cache_dir)
    version_dir <- file.path(cache_dir, as.character(row$version[[1L]]))
    dir.create(version_dir, recursive = TRUE, showWarnings = FALSE)
    destination <- file.path(version_dir, as.character(row$file_name[[1L]]))
    cached <- .spamtpdata_cached_file(row, cache_dir)
    if (!is.null(cached)) return(cached)

    url <- paste0(
        as.character(row$location_prefix[[1L]]),
        as.character(row$rdata_path[[1L]])
    )
    retries <- suppressWarnings(as.integer(retries)[1L])
    if (is.na(retries) || retries < 1L) retries <- 1L
    timeout <- suppressWarnings(as.numeric(timeout)[1L])
    if (!is.finite(timeout) || timeout < 1) timeout <- 1800
    old_timeout <- getOption("timeout")
    old_timeout_numeric <- suppressWarnings(as.numeric(old_timeout)[1L])
    if (!is.finite(old_timeout_numeric)) old_timeout_numeric <- 60
    options(timeout = max(old_timeout_numeric, timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    last_error <- NULL
    for (attempt in seq_len(retries)) {
        partial <- paste0(destination, ".part-", Sys.getpid())
        on.exit(unlink(partial), add = TRUE)
        result <- tryCatch(
            {
                download.file(
                    url,
                    destfile = partial,
                    method = "libcurl",
                    mode = "wb",
                    quiet = TRUE
                )
                if (!.spamtpdata_file_valid(partial, row)) {
                    stop("downloaded file failed its size or MD5 check")
                }
                if (file.exists(destination)) unlink(destination)
                if (!file.rename(partial, destination)) {
                    stop("could not move the verified file into the cache")
                }
                destination
            },
            error = function(error) {
                last_error <<- conditionMessage(error)
                unlink(partial)
                NULL
            }
        )
        if (!is.null(result)) return(result)
    }
    stop(
        "Failed to download verified SpaMTPData resource '", row$resource,
        "' after ", retries, " attempt(s): ", last_error,
        call. = FALSE
    )
}

#' List SpaMTP experiment resources
#'
#' @param version Data release version. `NULL` returns all versions.
#' @param category Optional resource category.
#'
#' @return A data frame describing the registered resources.
#' @export
#'
#' @examples
#' spaMTPDataResources()
spaMTPDataResources <- function(version = NULL, category = NULL) {
    manifest <- .spamtpdata_manifest()
    if (!is.null(version)) {
        version <- .spamtpdata_resolve_version(manifest, version)
        manifest <- manifest[manifest$version == version, , drop = FALSE]
    }
    if (!is.null(category)) {
        manifest <- manifest[manifest$category %in% category, , drop = FALSE]
    }
    rownames(manifest) <- NULL
    manifest
}

#' Retrieve a SpaMTP experiment resource
#'
#' Returns the representation actually stored in the published file; no
#' implicit conversion is performed. Release `1.1.0`, selected by default,
#' supplies seven native SpatialExperiment datasets. Use `version = "1.0.0"`
#' to request a historical release explicitly. Registry entries with
#' `requires_conversion = TRUE` are historical Seurat objects. Convert these
#' explicitly at the SpaMTP workflow boundary, or use [spaMTPExampleData()] for
#' a small native example that needs neither Seurat nor network access.
#' Local files and a verified download cache are checked before any Hub query.
#' Existing immutable files retain their original classes and checksums.
#'
#' @param resource Resource name; see [spaMTPDataResources()].
#' @param version Data release version or `"latest"`.
#' @param local_dir Optional directory containing downloaded source files.
#' @param hub Optional pre-created `ExperimentHub` object.
#' @param metadata Return only the registry row.
#' @param offline If `TRUE`, use local files or a verified cache only; never
#'   query ExperimentHub or the source URL.
#' @param fallback_url If `TRUE`, use the immutable source URL when the resource
#'   has not yet been ingested into ExperimentHub.
#' @param cache_dir Cache directory for source-URL downloads. Defaults to the
#'   platform-specific user cache returned by [tools::R_user_dir()].
#' @param timeout Download timeout in seconds for the source-URL fallback.
#' @param retries Number of verified download attempts.
#' @param verify Verify local files against the published size and checksum.
#'   Set to `FALSE` only for intentional development fixtures. Downloaded and
#'   cached files are always verified; loaded object classes are always checked.
#'
#' @return The requested experiment object or local file path. With
#'   `metadata = TRUE`, returns one registry row.
#' @export
#'
#' @examples
#' spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE)
spaMTPData <- function(resource, version = "latest", local_dir = NULL,
                       hub = NULL, metadata = FALSE, offline = FALSE,
                       fallback_url = TRUE, cache_dir = NULL, timeout = 1800,
                       retries = 3L, verify = TRUE) {
    if (length(resource) != 1L || is.na(resource) || !nzchar(trimws(resource))) {
        stop("resource must be one non-empty name.", call. = FALSE)
    }
    if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) {
        stop("verify must be TRUE or FALSE.", call. = FALSE)
    }
    manifest <- .spamtpdata_manifest()
    version <- .spamtpdata_resolve_version(manifest, version)
    key <- tolower(as.character(resource)[1L])
    rows <- manifest[
        tolower(manifest$resource) == key & manifest$version == version,
        , drop = FALSE
    ]
    if (nrow(rows) != 1L) {
        stop(
            "Unknown SpaMTPData resource '", resource, "' for version ",
            version, ". Use spaMTPDataResources() to list valid names.",
            call. = FALSE
        )
    }
    if (isTRUE(metadata)) return(rows)

    local_file <- .spamtpdata_local_file(rows, .spamtpdata_local_dir(local_dir))
    if (!is.null(local_file)) {
        .spamtpdata_verify_local(local_file, rows, verify)
        return(.spamtpdata_check_class(
            .spamtpdata_read_local(local_file, rows$dispatch_class), rows))
    }
    cached <- .spamtpdata_cached_file(rows, cache_dir)
    if (!is.null(cached)) {
        return(.spamtpdata_check_class(
            .spamtpdata_read_local(cached, rows$dispatch_class), rows))
    }
    if (isTRUE(offline)) {
        stop(
            "Resource '", resource, "' is not present in the configured local ",
            "directory or verified cache and offline = TRUE.",
            call. = FALSE
        )
    }

    hub_error <- NULL
    value <- tryCatch(
        {
            if (is.null(hub)) hub <- ExperimentHub()
            hits <- query(hub, c("SpaMTPData", rows$title))
            hit_metadata <- as.data.frame(mcols(hits))
            exact <- which(as.character(hit_metadata$title) == rows$title)
            if (!length(exact)) {
                stop("resource has not yet been ingested into ExperimentHub")
            }
            hits[[exact[[1L]]]]
        },
        error = function(error) {
            hub_error <<- conditionMessage(error)
            NULL
        }
    )
    if (!is.null(value)) return(.spamtpdata_check_class(value, rows))
    if (isTRUE(fallback_url)) {
        path <- .spamtpdata_download(
            rows,
            cache_dir = cache_dir,
            timeout = timeout,
            retries = retries
        )
        return(.spamtpdata_check_class(
            .spamtpdata_read_local(path, rows$dispatch_class), rows))
    }
    stop(
        "ExperimentHub could not provide '", rows$title, "': ", hub_error,
        ". Configure a local resource directory or set fallback_url = TRUE.",
        call. = FALSE
    )
}

#' Alias for explicit SpaMTPData resource retrieval
#'
#' @inheritParams spaMTPData
#' @return The value returned by [spaMTPData()].
#' @export
#'
#' @examples
#' spaMTPDataResource("mouse_brain_dhb_striatum", metadata = TRUE)
spaMTPDataResource <- function(resource, version = "latest", local_dir = NULL,
                               hub = NULL, metadata = FALSE, offline = FALSE,
                               fallback_url = TRUE, cache_dir = NULL,
                               timeout = 1800, retries = 3L, verify = TRUE) {
    spaMTPData(
        resource = resource,
        version = version,
        local_dir = local_dir,
        hub = hub,
        metadata = metadata,
        offline = offline,
        fallback_url = fallback_url,
        cache_dir = cache_dir,
        timeout = timeout,
        retries = retries,
        verify = verify
    )
}

#' Report the available SpaMTPData releases
#'
#' @return A character vector of data releases, newest first.
#' @export
#'
#' @examples
#' spaMTPDataVersion()
spaMTPDataVersion <- function() {
    versions <- unique(as.character(.spamtpdata_manifest()$version))
    rev(as.character(sort(package_version(versions))))
}
