test_that("the native example is reproducible and uses public S4 containers", {
    before <- intersect(c("Seurat", "SeuratObject"), loadedNamespaces())
    object <- spaMTPExampleData()
    expect_s4_class(object, "SpatialExperiment")
    expect_true(methods::validObject(object))
    expect_identical(object, spaMTPExampleData())
    expect_equal(dim(object), c(6L, 12L))
    expect_identical(SummarizedExperiment::assayNames(object), "counts")
    expect_true(all(is.finite(SummarizedExperiment::rowData(object)$mz)))
    expect_true(all(SummarizedExperiment::assay(object) >= 0))
    expect_named(as.data.frame(SpatialExperiment::spatialCoords(object)), c("x", "y"))
    expect_false(any(c("x", "y", "x_coord", "y_coord") %in%
        colnames(SummarizedExperiment::colData(object))))
    expect_true(S4Vectors::metadata(object)$example$synthetic)
    expect_identical(intersect(c("Seurat", "SeuratObject"), loadedNamespaces()), before)
})

test_that("alternative experiments follow pixel subsetting", {
    object <- spaMTPExampleData()
    paired <- SingleCellExperiment::altExp(object, "transcriptome")
    expect_equal(dim(paired), c(5L, 12L))
    expect_identical(colnames(paired), colnames(object))
    selected <- object[, object$region == "core"]
    expect_identical(colnames(SingleCellExperiment::altExp(selected, "transcriptome")),
        colnames(selected))
    expect_equal(nrow(SpatialExperiment::spatialCoords(selected)), 6L)
    expect_length(SingleCellExperiment::altExpNames(spaMTPExampleData(FALSE)), 0L)
    expect_error(spaMTPExampleData(NA), "TRUE or FALSE")
})
