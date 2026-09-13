# One-time native resource preparation

`prepare-native-resources.R` turns the seven published Seurat archives into
standalone SpatialExperiment RDS files. Conversion is a **maintainer step**;
neither SpaMTPData nor later analysis of these RDS files needs Seurat.
The large files belong in an external data release, not Git or the R source
package. Resource release `1.1.0` is published in
[Zenodo record 22733262](https://zenodo.org/records/22733262) and is the default
in SpaMTPData >= 0.99.4. Running the preparation recipe does not automatically
upload files or register them with ExperimentHub.

## Prepare

Install the coordinated SpaMTP >= 0.99.4 software and its dependencies, plus
Seurat and SeuratObject for reading historical image classes. The conversion
uses public accessors and the SpaMTP adapters. It does not access slots directly.

Run from a SpaMTPData source checkout:

```sh
Rscript inst/scripts/prepare-native-resources.R \
  /path/to/official-archives /path/to/native-staging 1.1.0
```

The input directory may contain a `1.0.0` subdirectory. Alternatively, the first
argument may be a CSV with `resource`, `version` and absolute `path` columns,
allowing files already present in different locations to be reused. All input
paths are checked against the **package's published registry**, not checksums
supplied by the input CSV. Modified files are rejected. The recipe never calls
an upload API and does not update published resource or Hub metadata.

The conversion plan is `inst/manifest/native_conversion_plan.csv`. Primary
assays and images are explicit, including the `image` view for the human brain
and the supplied `centroids` view for simulated Xenium. No pairing is inferred
between independently measured Xenium/MALDI or Visium/MSI datasets.

Each output has a `_spe.rds` suffix. The accompanying
`native_resource_manifest.csv` records dimensions, layers, paired assays, source
file identifiers/checksums and native-file checksums. Metadata in each RDS also
records the preparation code hashes, software versions and omitted components.
An existing file with different provenance is not overwritten.

## Retained and omitted data

- Exact sparse `counts` matrices, plus complete `data` layers when present.
  `data` is not renamed to `logcounts`: its numerical values and interpretation
  remain those of the source archive. Scaling/PCA can be computed natively.
- All feature metadata, plus numeric m/z for metabolomics. Synthetic m/z and
  gene names in the simulated datasets remain explicitly labelled synthetic.
- Pixel metadata and identities, including human-brain protein intensities and
  non-syntactic metabolite column names. Shared pixel metadata is stored once.
- Explicitly selected unscaled image-centroid coordinates in `spatialCoords()`.
  Redundant x/y metadata is removed; original source files remain available.
- Paired `SPT` matrices in `altExp()` for DHB striatum and the human brain,
  with exactly matched pixel names. Portable stored annotation tables remain
  available, including legacy `db_3` with its historical annotation semantics.

Optical rasters, segmentation polygons, molecule coordinates, `scale.data`,
graphs, reductions, Seurat command history and Seurat-specific weight objects
are omitted. Native data therefore support expression/annotation/spatial
workflows, but do not reproduce historical polygon-overlap mapping or WNN
embeddings. Add an optical image or explicit sf polygons separately if needed.

Original sample labels are retained. Do not assume two independently imported
experiments have matching `sample_id` values or coordinate units; confirm and
align them before spatial mapping.

## Use the resulting RDS directly

```r
spe <- readRDS("/path/to/native-staging/1.1.0/mouse_brain_dhb_striatum_spe.rds")
SummarizedExperiment::assayNames(spe)
SpatialExperiment::spatialCoords(spe)
SingleCellExperiment::altExp(spe, "SPT")
spe <- SpaMTP::normalizeSMData(spe, verbose = FALSE)
```

No Seurat conversion is needed at this point. To retrieve the published native
file instead, use `SpaMTPData::spaMTPData("mouse_brain_dhb_striatum",
version = "1.1.0")`. Explicit `version = "1.0.0"` requests retain the historical
archives. Do not place a native file under an old archive's filename and
checksum. Publishing files and switching the registry remain separate release
steps; publication does not itself complete ExperimentHub ingestion.

Publication is coordinated by SpaMTPdb's `prepare-shared-release.R` and
`upload-shared-release.R` recipes in the existing Zenodo version family. The
experimental RDS can share storage with annotation databases without becoming
members of SpaMTPdb's runtime registry. `register-native-release.R` in this
directory then verifies a published release and fresh native reads, and generates
candidate SpaMTPData registry/Hub tables containing all 18 resources per release.
The seven native replacements do not make the eleven unchanged resources
disappear from the default version. No separate SpaMTPData Zenodo family is needed.

Coordinate-accessor references:
[SeuratObject tissue coordinates](https://satijalab.github.io/seurat-object/reference/GetTissueCoordinates.html),
[centroid access](https://satijalab.github.io/seurat-object/reference/Centroids-methods.html).
