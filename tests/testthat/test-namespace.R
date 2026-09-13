test_that("every Imports dependency has explicit namespace imports", {
    declared <- utils::packageDescription("SpaMTPData", fields = "Imports")
    dependencies <- trimws(gsub("\\s*\\([^)]*\\)", "",
        strsplit(declared, ",")[[1L]]))
    imports <- getNamespaceImports("SpaMTPData")

    expect_true(all(dependencies %in% names(imports)))
    for (dependency in dependencies) {
        symbols <- imports[[dependency]]
        expect_type(symbols, "character")
        expect_gt(length(symbols), 0L)
    }
})

test_that("S4 replacement accessors are explicitly imported", {
    imports <- getNamespaceImports("SpaMTPData")
    expect_true(all(c("metadata", "metadata<-") %in% imports$S4Vectors))
    expect_true("altExp<-" %in% imports$SingleCellExperiment)
})
