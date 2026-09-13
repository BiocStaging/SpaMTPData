test_that("native resource provenance includes declared species without modifying assays", {
    object <- spaMTPExampleData()
    before <- S4Vectors::metadata(object)
    for (resource in c("human_brain_integrated", "mouse_brain_dhb_striatum")) {
        row <- spaMTPData(resource, metadata = TRUE)
        loaded <- .spamtpdata_check_class(object, row)
        expected <- if (resource == "human_brain_integrated") "Homo sapiens" else "Mus musculus"
        expect_identical(S4Vectors::metadata(loaded)$SpaMTPData$organism, expected)
        expect_identical(S4Vectors::metadata(loaded)$SpaMTPData$md5, row$md5)
        expect_identical(S4Vectors::metadata(loaded)$example, before$example)
        expect_identical(SummarizedExperiment::assay(loaded), SummarizedExperiment::assay(object))
        paired <- SingleCellExperiment::altExp(loaded, "transcriptome")
        expect_identical(S4Vectors::metadata(paired)$SpaMTPData,
                         S4Vectors::metadata(loaded)$SpaMTPData)
        S4Vectors::metadata(paired)$SpaMTPData <- NULL
        expect_identical(paired, SingleCellExperiment::altExp(object, "transcriptome"))
        expect_identical(SpatialExperiment::spatialCoords(loaded), SpatialExperiment::spatialCoords(object))
    }
})

test_that("species provenance survives nested and unnamed altExp extraction", {
    object <- spaMTPExampleData()
    child <- SingleCellExperiment::altExp(object, "transcriptome")
    nested <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(child)))
    S4Vectors::metadata(nested)$note <- "retain me"
    SingleCellExperiment::altExp(child, "nested") <- nested
    SingleCellExperiment::altExp(object, "transcriptome") <- child
    SingleCellExperiment::altExp(object, "second") <- nested
    expect_warning(SingleCellExperiment::altExpNames(object) <- c("", ""), "empty strings")
    original_names <- SingleCellExperiment::altExpNames(object)
    row <- spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE)
    loaded <- .spamtpdata_check_class(object, row)
    first <- SingleCellExperiment::altExp(loaded, 1L)
    deepest <- SingleCellExperiment::altExp(first, "nested")
    second <- SingleCellExperiment::altExp(loaded, 2L)
    for (x in list(first, deepest, second)) {
        provenance <- S4Vectors::metadata(x)$SpaMTPData
        expect_identical(provenance$organism, "Mus musculus")
        expect_identical(provenance$taxonomy_id, 10090L)
        expect_identical(provenance$genome, "GRCm39")
        expect_identical(provenance$md5, row$md5)
    }
    expect_identical(S4Vectors::metadata(deepest)$note, "retain me")
    expect_identical(SingleCellExperiment::altExpNames(loaded), original_names)
    S4Vectors::metadata(deepest)$SpaMTPData <- NULL
    expect_identical(deepest, nested)
})

test_that("resource species does not relabel explicitly declared child experiments", {
    object <- spaMTPExampleData()
    child <- SingleCellExperiment::altExp(object, "transcriptome")
    S4Vectors::metadata(child)$organism <- "Homo sapiens"
    S4Vectors::metadata(child)$genome <- "GRCh38"
    SingleCellExperiment::altExp(object, "transcriptome") <- child
    loaded <- .spamtpdata_check_class(object,
        spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE))
    provenance <- S4Vectors::metadata(SingleCellExperiment::altExp(loaded, "transcriptome"))$SpaMTPData
    expect_identical(provenance$organism, "Homo sapiens")
    expect_identical(provenance$genome, "GRCh38")
    expect_null(provenance$taxonomy_id)
    expect_identical(S4Vectors::metadata(loaded)$SpaMTPData$organism, "Mus musculus")
})
