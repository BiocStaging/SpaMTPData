# Maintainer-only helpers: Seurat is used once to produce portable native RDS.
# Runtime SpaMTPData resource retrieval never sources this file.

nativePortableMetadata <- function(value) {
    if (inherits(value, "data.frame")) {
        value <- as.data.frame(value)
        value[] <- lapply(value, nativePortableMetadata)
        return(value)
    }
    if (isS4(value) || is.environment(value) || is.function(value) ||
        typeof(value) %in% c("externalptr", "weakref")) {
        stop("Non-portable value in stored analysis metadata: ",
             paste(class(value), collapse = "/"), call. = FALSE)
    }
    if (is.list(value)) return(lapply(value, nativePortableMetadata))
    if (is.atomic(value) || is.null(value)) return(value)
    stop("Unsupported stored metadata type: ", typeof(value), call. = FALSE)
}

nativeSourceCoordinates <- function(x, image) {
    if (!image %in% SeuratObject::Images(x)) {
        stop("The conversion plan names an unavailable image: ", image)
    }
    source <- x[[image]]
    boundary <- NA_character_
    if (inherits(source, "FOV") &&
        "centroids" %in% SeuratObject::Boundaries(source)) {
        boundary <- "centroids"
        values <- SeuratObject::GetTissueCoordinates(
            source[[boundary]], full = FALSE)
    } else if (inherits(source, "VisiumV1")) {
        values <- SeuratObject::GetTissueCoordinates(source, scale = NULL)
    } else {
        values <- SeuratObject::GetTissueCoordinates(source, full = FALSE)
    }
    axes <- c(x = "x", y = "y")
    if (!all(axes %in% names(values)) &&
        all(c("imagecol", "imagerow") %in% names(values))) {
        axes <- c(x = "imagecol", y = "imagerow")
        values$x <- values$imagecol
        values$y <- values$imagerow
    }
    if (!all(c("x", "y") %in% names(values))) {
        stop("Image coordinates need explicit x/y axis mapping: ", image)
    }
    if (!"cell" %in% names(values)) {
        values$cell <- if ("ID" %in% names(values)) values$ID else rownames(values)
    }
    list(values = values[, c("x", "y", "cell")], image = image,
         image_class = class(source), boundary = boundary,
         axes = axes, units = "original unscaled source-image coordinates")
}

nativeCopyDataLayer <- function(target, source, assay) {
    available <- SeuratObject::Layers(source[[assay]], search = NA)
    if (!"counts" %in% available || any(grepl("^(counts|data)\\.", available))) {
        stop("The recipe requires complete, unsplit counts/data layers: ", assay)
    }
    if ("data" %in% available) {
        value <- SeuratObject::LayerData(source, assay = assay, layer = "data")
        if (!setequal(rownames(value), rownames(target)) ||
            !setequal(colnames(value), colnames(target))) {
            stop("The data layer does not cover the complete counts matrix: ", assay)
        }
        SummarizedExperiment::assay(target, "data") <-
            value[rownames(target), colnames(target), drop = FALSE]
    }
    target
}

nativeMassMetadata <- function(object) {
    features <- SummarizedExperiment::rowData(object)
    masses <- NULL
    for (column in c("mz", "raw_mz")) {
        if (column %in% colnames(features)) {
            candidate <- suppressWarnings(as.numeric(as.character(features[[column]])))
            if (length(candidate) == nrow(object) && all(is.finite(candidate))) {
                masses <- candidate
                break
            }
        }
    }
    if (is.null(masses) && all(grepl("^mz-", rownames(object)))) {
        masses <- suppressWarnings(as.numeric(sub("^mz-", "", rownames(object))))
    }
    if (is.null(masses) || length(masses) != nrow(object) ||
        any(!is.finite(masses)) || any(masses <= 0)) {
        stop("Metabolite features require unambiguous positive numeric m/z values.")
    }
    features$mz <- masses
    SummarizedExperiment::rowData(object) <- features
    object
}

prepareNativeResource <- function(x, specification) {
    assay <- specification$primary_assay
    alternatives <- strsplit(specification$alternative_assays, ";", fixed = TRUE)[[1L]]
    alternatives <- alternatives[nzchar(alternatives)]
    coordinates <- nativeSourceCoordinates(x, specification$image)
    object <- SpaMTP::seuratToSpatialExperiment(x, assay = assay,
        layer = "counts", includeAltExps = FALSE, coordinates = coordinates$values)
    object <- nativeCopyDataLayer(object, x, assay)
    if (specification$modality == "metabolome") object <- nativeMassMetadata(object)
    retained <- S4Vectors::metadata(object)
    retained <- retained[setdiff(names(retained), "polygon_coordinates")]
    retained <- nativePortableMetadata(retained)
    S4Vectors::metadata(object) <- retained

    # Pixel metadata is shared once in the primary experiment. Authoritative
    # coordinates live in spatialCoords(), not stale Seurat metadata columns.
    pixels <- SummarizedExperiment::colData(object)
    removedCoordinates <- intersect(c("x", "y", "x_coord", "y_coord",
        "x_centroid", "y_centroid"), colnames(pixels))
    SummarizedExperiment::colData(object) <-
        pixels[, setdiff(colnames(pixels), removedCoordinates), drop = FALSE]
    for (alternative in alternatives) {
        experiment <- SpaMTP::seuratToSingleCellExperiment(x, assay = alternative,
            layer = "counts", includeAltExps = FALSE)
        if (!setequal(colnames(experiment), colnames(object))) {
            stop("Alternative assay is not paired at the same pixels: ", alternative)
        }
        experiment <- nativeCopyDataLayer(experiment, x, alternative)
        experiment <- experiment[, colnames(object), drop = FALSE]
        SummarizedExperiment::colData(experiment) <-
            S4Vectors::DataFrame(row.names = colnames(experiment))
        S4Vectors::metadata(experiment) <- list(source_assay = alternative)
        SingleCellExperiment::altExp(object, alternative) <- experiment
    }
    sourceAssays <- SeuratObject::Assays(x)
    sourceLayers <- lapply(c(assay, alternatives), function(name) {
        SeuratObject::Layers(x[[name]], search = NA)
    })
    names(sourceLayers) <- c(assay, alternatives)
    S4Vectors::metadata(object)$native_conversion <- list(
        recipe_version = "1", resource = specification$resource,
        modality = specification$modality, synthetic = specification$synthetic,
        primary_assay = assay, alternative_assays = alternatives,
        layers = c("counts", if ("data" %in% SummarizedExperiment::assayNames(object)) "data"),
        source_layers = sourceLayers,
        coordinates = coordinates[setdiff(names(coordinates), "values")],
        omitted_assays = setdiff(sourceAssays, c(assay, alternatives)),
        omitted_layers = lapply(sourceLayers, setdiff, c("counts", "data")),
        omitted_components = c("optical rasters", "segmentation polygons",
            "molecule coordinates", "graphs", "reductions", "Seurat commands"),
        omitted_metadata = setdiff(union(names(SeuratObject::Misc(x)),
            SeuratObject::Tool(x)), names(retained)),
        omitted_coordinate_columns = removedCoordinates,
        data_layer = "preserved verbatim; not relabelled or assumed log-normalized")
    methods::validObject(object)
    object
}

validateNativeResource <- function(object, source, specification) {
    stopifnot(methods::is(object, "SpatialExperiment"), methods::validObject(object),
        length(SingleCellExperiment::reducedDimNames(object)) == 0L,
        nrow(SpatialExperiment::imgData(object)) == 0L)
    names <- c(specification$primary_assay,
        SingleCellExperiment::altExpNames(object))
    for (name in names) {
        experiment <- if (name == specification$primary_assay) object else
            SingleCellExperiment::altExp(object, name)
        for (layer in SummarizedExperiment::assayNames(experiment)) {
            original <- SeuratObject::LayerData(source, assay = name, layer = layer)
            original <- original[rownames(experiment), colnames(experiment), drop = FALSE]
            stopifnot(identical(SummarizedExperiment::assay(experiment, layer), original))
        }
        features <- source[[name]][[]]
        converted <- as.data.frame(SummarizedExperiment::rowData(experiment), optional = TRUE)
        # Empty feature metadata can differ only in NULL versus character(0)
        # column-name attributes after DataFrame conversion. Counts supply the
        # feature IDs; there are no metadata column values to compare.
        if (ncol(features)) {
            stopifnot(identical(converted[, names(features), drop = FALSE],
                as.data.frame(features[rownames(experiment), , drop = FALSE])))
        }
    }
    originalPixels <- source[[]]
    for (name in setdiff(names(originalPixels), c("x", "y", "x_coord", "y_coord",
        "x_centroid", "y_centroid"))) {
        stopifnot(identical(SummarizedExperiment::colData(object)[[name]],
            originalPixels[colnames(object), name]))
    }
    coordinates <- nativeSourceCoordinates(source, specification$image)$values
    index <- match(colnames(object), as.character(coordinates$cell))
    stopifnot(!anyNA(index), !anyDuplicated(coordinates$cell),
        identical(unname(SpatialExperiment::spatialCoords(object)),
            unname(as.matrix(coordinates[index, c("x", "y")]))))
    invisible(TRUE)
}
