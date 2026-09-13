resourceFixture <- function(directory, value = list(dataset = "synthetic")) {
    path <- file.path(directory, "fixture.rds")
    saveRDS(value, path)
    data.frame(resource = "fixture", version = "1.0.0", category = "test",
        title = "SpaMTPData_fixture_1.0.0", file_name = basename(path),
        location_prefix = paste0("file://", normalizePath(directory), "/"),
        rdata_path = basename(path), r_data_class = paste(class(value), collapse = "/"),
        dispatch_class = "Rds", bytes = file.info(path)$size,
        md5 = unname(tools::md5sum(path)), stringsAsFactors = FALSE)
}
