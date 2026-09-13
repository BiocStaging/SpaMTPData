#' Construct a small native spatial multi-omics example
#'
#' A deterministic, synthetic example for tutorials and interoperability tests.
#' This is not a biological dataset or a downloaded ExperimentHub resource.
#' The primary experiment contains six m/z features measured at twelve pixels;
#' an optional paired transcriptome has five genes at the same pixels.
#' All data use existing Bioconductor containers and public accessors. Neither
#' SpaMTP nor Seurat is required to construct or inspect the example.
#'
#' @param paired Include a transcriptome in `altExp(x, "transcriptome")`.
#'
#' @return A `SpatialExperiment` with intensities in `assay(x, "counts")`, m/z
#'   values in `rowData(x)`, regions in `colData(x)`, pixel x/y coordinates in
#'   `spatialCoords(x)` and synthetic-data provenance in `metadata(x)`.
#' @export
#' @examples
#' spe <- spaMTPExampleData()
#' SummarizedExperiment::assayNames(spe)
#' SummarizedExperiment::rowData(spe)$mz
#' SpatialExperiment::spatialCoords(spe)
#' SingleCellExperiment::altExpNames(spe)
spaMTPExampleData <- function(paired = TRUE) {
    if (!is.logical(paired) || length(paired) != 1L || is.na(paired)) {
        stop("paired must be TRUE or FALSE.", call. = FALSE)
    }
    pixels <- paste0("pixel", seq_len(12L))
    counts <- outer(seq_len(6L), seq_len(12L), function(feature, pixel) {
        5 + (feature * pixel + pixel^2 + 3 * feature^2) %% 23 +
            20 * (pixel <= 6 & feature <= 3)
    })
    dimnames(counts) <- list(paste0("peak", seq_len(6L)), pixels)
    # Illustrative protonated/sodiated masses plus one unmatched feature.
    mz <- c(149.0444498, 171.0263941, 181.0706646, 203.0526088,
            166.0862550, 500.05)
    object <- SpatialExperiment(
        assays = list(counts = counts),
        rowData = DataFrame(mz = mz, raw_mz = mz),
        colData = DataFrame(
            region = factor(rep(c("edge", "core"), each = 6L)),
            row.names = pixels),
        spatialCoords = cbind(x = rep(0:3, each = 3L), y = rep(0:2, 4L)),
        sample_id = "synthetic-section1"
    )
    if (paired) {
        transcriptome <- outer(seq_len(5L), seq_len(12L), function(gene, pixel) {
            1 + (11 * gene + pixel^2 + gene * pixel) %% 19
        })
        dimnames(transcriptome) <- list(paste0("gene", seq_len(5L)), pixels)
        altExp(object, "transcriptome") <-
            SingleCellExperiment(
                assays = list(counts = transcriptome))
    }
    metadata(object)$example <- list(
        synthetic = TRUE, generator = "SpaMTPData::spaMTPExampleData",
        generator_version = "1", coordinate_units = "arbitrary pixel grid")
    object
}
