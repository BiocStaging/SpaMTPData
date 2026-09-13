# SpaMTPData

SpaMTPData >= 0.99.5 attaches declared species and resource provenance to loaded
native experiments and their nested `altExp()` objects in
`metadata(x)$SpaMTPData`. Species checks therefore survive extracting a paired
transcriptome, while explicit child species declarations are retained.
SpaMTP >= 0.99.5 can use these
fields when annotating genes with `annotateGeneIdentifiers()`, whose HGNC
reference is provided by SpaMTPdb >= 0.99.4. Human HGNC mapping is rejected for
experiments declared as mouse; it does not perform orthology conversion.
The published experiment files, their hashes and resource versions are preserved.

`SpaMTPData` provides ExperimentHub access to the datasets used by
[`SpaMTP`](https://github.com/SpaMTP-project/SpaMTP) tutorials and
tests. SpatialExperiment datasets, Cardinal experiments and auxiliary files
remain outside Git and are downloaded only when requested. ExperimentHub is preferred; before Hub
ingestion is complete, the same API falls back to immutable Zenodo URLs with a
30-minute timeout, retries, and size/MD5 validation.

```r
# Native, deterministic synthetic example: no download or Seurat dependency.
spe <- SpaMTPData::spaMTPExampleData()
SummarizedExperiment::rowData(spe)$mz
SingleCellExperiment::altExpNames(spe)

# Inspect a published resource before downloading its original representation.
SpaMTPData::spaMTPData("mouse_brain_dhb_striatum", metadata = TRUE)
```

The coordinated API uses `spaMTPData()`, `spaMTPDataResource()`,
`spaMTPDataResources()` and `spaMTPDataVersion()`. This replaces the historical
PascalCase exports. Use SpaMTPData >= 0.99.4 with SpaMTP >= 0.99.4 and
SpaMTPdb >= 0.99.3 for the coordinated release. Package versions are distinct
from resource releases `1.1.0` and `1.0.0`.

The default resource release, `1.1.0`, provides seven native SpatialExperiment
datasets from [Zenodo record 22733262](https://zenodo.org/records/22733262).
The historical `1.0.0` release remains available explicitly; its seven Seurat
entries retain `requires_conversion = TRUE`, their original URLs and hashes.
Historical objects are not silently converted or relabelled.

The [one-time native preparation recipe](inst/scripts/README-native-resources.md)
extracts the required data from all seven archives into standalone
SpatialExperiment RDS. The published `1.1.0` files retain counts/data, annotations,
source-image centroids and paired transcriptomes, and can be used directly:

```r
spe <- SpaMTPData::spaMTPData("mouse_brain_dhb_striatum", version = "1.1.0")
SingleCellExperiment::altExpNames(spe) # SPT
```

Preparation needs optional Seurat support once; reading these native RDS does
not. The prepared files have their own provenance and checksum manifest; they
are not uploaded by the recipe or silently substituted for old published URLs.

For offline development, set `options(SpaMTPData.resource_dir = "...")` to a
directory containing files named as listed in the resource registry. A
version-specific subdirectory is preferred. Local files are checked against
published sizes/checksums by default; `verify = FALSE` is for deliberate
development fixtures only, and never disables class validation.
Verified source-download caches also work with `offline = TRUE`, without a Hub
connection. Cache entries are versioned to avoid filename collisions.

Each release contains 18 resources: the current native snapshot and the
historical snapshot together have 36 registry entries, covering mouse brain, mouse bladder, import examples,
simulated single-cell multi-omics, annotation refinement, a large-data ROI, and
the human brain demonstration.

`SpaMTPData` owns experiment resources, `SpaMTPdb` owns annotation tables and
SpaMTP owns analysis/conversion. Neither data package imports SpaMTP or Seurat.
Physical files share the existing SpaMTPdb Zenodo version family; this does
not require a separate SpaMTPData deposition family or a dependency between
the two R packages. The shared publication recipes live in SpaMTPdb >= 0.99.3.
The native synthetic example is separate from the published-resource registry;
it is a teaching fixture, not an additional biological dataset.

After a shared native release is publicly available, run the registration recipe
from a SpaMTPData >= 0.99.3 checkout:

```sh
Rscript inst/scripts/register-native-release.R \
  /path/to/shared-staging PUBLISHED_RECORD_ID /path/to/new-registration-output
```

This rejects draft IDs and other Zenodo version families, verifies file sizes
and MD5 against the prepared native manifest, and reads all seven public files
in a fresh temporary cache as SpatialExperiment without Seurat. It writes
candidate `resource_manifest.csv` and `metadata.csv` files to the specified
output directory; it does not modify this package's live registry or ingest Hub
records. Review and validate the candidate tables before installing them in
`inst/manifest` and `inst/extdata`, respectively.

The new experiment release is a complete snapshot: seven native replacements
plus the eleven unchanged resources at their existing immutable locations.
Historical 1.0.0 entries remain available. Publishing a new resource does not
silently replace the source class or checksum of an older release.

If you use these resources, cite the SpaMTP paper:
[Causer, Lu, Kriel *et al.*, Nature Methods (2026)](https://doi.org/10.1038/s41592-026-03140-8).
