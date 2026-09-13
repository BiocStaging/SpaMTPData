test_that("resource registry is versioned and categorized", {
    resources <- spaMTPDataResources()
    expect_true(all(c("resource", "version", "category") %in% names(resources)))
    expect_true("mouse_brain_dhb_striatum" %in% resources$resource)
    expect_true(all(resources$bytes > 0))
    expect_true(all(grepl("^[0-9a-f]{32}$", resources$md5)))
    expect_identical(spaMTPDataVersion(), c("1.1.0", "1.0.0"))
    expect_equal(resources$requires_conversion, resources$r_data_class == "Seurat")
})

test_that("the native release is the complete default snapshot with preserved history", {
    latest <- spaMTPDataResources(version = "latest")
    old <- spaMTPDataResources(version = "1.0.0")
    expect_identical(unique(latest$version), "1.1.0")
    expect_equal(nrow(latest), 18L)
    expect_equal(nrow(old), 18L)
    expect_setequal(latest$resource, old$resource)
    expect_false(any(latest$requires_conversion))
    native <- latest$r_data_class == "SpatialExperiment"
    expect_equal(sum(native), 7L)
    expect_true(all(latest$location_prefix[native] ==
        "https://zenodo.org/api/records/22733262/files/"))
    expect_true(all(latest$source_url[native] ==
        "https://zenodo.org/records/22733262"))
    old <- old[match(latest$resource, old$resource), , drop = FALSE]
    expect_true(all(old$requires_conversion[native]))
    expect_true(all(latest$md5[native] != old$md5[native]))
    for (column in c("file_name", "location_prefix", "rdata_path", "source_url",
        "r_data_class", "dispatch_class", "bytes", "md5")) {
        expect_identical(latest[[column]][!native], old[[column]][!native])
    }
    metadata <- spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE)
    expect_identical(metadata$version, "1.1.0")
    expect_identical(metadata$r_data_class, "SpatialExperiment")
    historical <- spaMTPData("mouse_brain_dhb_striatum",
        version = "1.0.0", metadata = TRUE)
    expect_identical(historical$r_data_class, "Seurat")
    expect_identical(historical$md5, "8fca750fc14909e638450fb5c5a0512b")
})

test_that("Hub metadata matches both releases and retains original source attribution", {
    registry <- spaMTPDataResources()
    metadata <- read.csv(system.file("extdata", "metadata.csv",
        package = "SpaMTPData"), stringsAsFactors = FALSE)
    index <- match(registry$title, metadata$Title)
    expect_equal(nrow(metadata), 36L)
    expect_false(anyNA(index))
    expect_false(anyDuplicated(metadata$Title) > 0L)
    for (pair in list(c("r_data_class", "RDataClass"),
        c("dispatch_class", "DispatchClass"), c("rdata_path", "RDataPath"),
        c("location_prefix", "Location_Prefix"))) {
        expect_identical(registry[[pair[1L]]], metadata[[pair[2L]]][index])
    }
    native <- registry$version == "1.1.0" &
        registry$r_data_class == "SpatialExperiment"
    source <- spaMTPDataResources(version = "1.0.0")
    source <- source[match(registry$resource[native], source$resource), ]
    expect_identical(metadata$SourceUrl[index[native]], source$source_url)
    expect_identical(metadata$SourceVersion[index[native]], source$version)
})

test_that("local RDS resources load without a Hub connection", {
    path <- tempfile("spamtpdata-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    fixture <- list(dataset = "mouse brain")
    row <- resourceFixture(path, fixture)
    local_mocked_bindings(.spamtpdata_manifest = function() row)
    observed <- spaMTPData(
        "fixture",
        local_dir = path,
        offline = TRUE
    )
    expect_identical(observed, fixture)
})

test_that("metadata lookup never downloads a resource", {
    metadata <- spaMTPData("simulated_xenium", metadata = TRUE)
    expect_identical(metadata$category, "simulated_multiomics")
    expect_identical(metadata$dispatch_class, "Rds")
})

test_that("source fallback verifies and reuses its cache", {
    source_dir <- tempfile("spamtpdata-source-")
    cache_dir <- tempfile("spamtpdata-cache-")
    dir.create(source_dir)
    dir.create(cache_dir)
    on.exit(unlink(c(source_dir, cache_dir), recursive = TRUE), add = TRUE)

    source_file <- file.path(source_dir, "fixture.rds")
    saveRDS(list(value = 42), source_file)
    row <- data.frame(
        resource = "fixture",
        version = "1.0.0",
        file_name = "fixture.rds",
        location_prefix = paste0("file://", normalizePath(source_dir), "/"),
        rdata_path = "fixture.rds",
        bytes = file.info(source_file)$size,
        md5 = unname(tools::md5sum(source_file)),
        stringsAsFactors = FALSE
    )

    first <- .spamtpdata_download(row, cache_dir = cache_dir, retries = 1L)
    second <- .spamtpdata_download(row, cache_dir = cache_dir, retries = 1L)
    expect_identical(first, second)
    expect_true(.spamtpdata_file_valid(second, row))
})

test_that("offline access reuses a verified versioned cache without a local file", {
    path <- tempfile("spamtpdata-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    fixture <- list(value = 42)
    row <- resourceFixture(path, fixture)
    cache <- file.path(path, "cache")
    .spamtpdata_download(row, cache_dir = cache, retries = 1L)
    local_mocked_bindings(.spamtpdata_manifest = function() row)
    expect_identical(spaMTPData("fixture", local_dir = file.path(path, "empty"),
        cache_dir = cache, offline = TRUE), fixture)
    expect_identical(spaMTPDataResource("fixture", local_dir = file.path(path, "empty"),
        cache_dir = cache, offline = TRUE), fixture)
    saveRDS("corrupt", file.path(cache, "1.0.0", "fixture.rds"))
    expect_error(spaMTPData("fixture", local_dir = file.path(path, "empty"),
        cache_dir = cache, offline = TRUE), "offline = TRUE")
})

test_that("local checksums and classes are enforced independently", {
    path <- tempfile("spamtpdata-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    row <- resourceFixture(path)
    local_mocked_bindings(.spamtpdata_manifest = function() row)
    saveRDS(list(changed = TRUE), file.path(path, row$file_name))
    expect_error(spaMTPData("fixture", local_dir = path, offline = TRUE), "MD5")
    expect_identical(spaMTPData("fixture", local_dir = path, offline = TRUE,
        verify = FALSE), list(changed = TRUE))
    saveRDS(42, file.path(path, row$file_name))
    expect_error(spaMTPData("fixture", local_dir = path, offline = TRUE,
        verify = FALSE), "expected list")
    expect_error(spaMTPData(character()), "one non-empty name")
    expect_error(spaMTPData("fixture", verify = NA), "TRUE or FALSE")
})

test_that("different versions with one filename do not overwrite each other", {
    path <- tempfile("spamtpdata-")
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    row1 <- resourceFixture(path, list(version = 1))
    cache <- file.path(path, "cache")
    first <- .spamtpdata_download(row1, cache, retries = 1L)
    row2 <- resourceFixture(path, list(version = 2))
    row2$version <- "2.0.0"
    second <- .spamtpdata_download(row2, cache, retries = 1L)
    expect_false(identical(first, second))
    expect_true(.spamtpdata_file_valid(first, row1))
    expect_true(.spamtpdata_file_valid(second, row2))
    expect_identical(readRDS(first)$version, 1)
    expect_identical(readRDS(second)$version, 2)
})

test_that("exports are camelCase and only Hub metadata is in extdata", {
    exports <- getNamespaceExports("SpaMTPData")
    expect_true(all(grepl("^[a-z][A-Za-z0-9]*$", exports)))
    expect_identical(list.files(system.file("extdata", package = "SpaMTPData"),
        pattern = "\\.csv$"), "metadata.csv")
})
