# SpaMTPData 0.99.5

* Loaded native experiments now retain resource version, checksum, organism,
  taxonomy ID and genome in `metadata(x)$SpaMTPData`. Analysis workflows can
  validate gene-reference species without guessing from gene capitalization.
  Published resource files, versions, assays and coordinates are unchanged.
* Propagated this provenance to nested alternative experiments using public
  accessors, preserving species checks after `altExp()` extraction. Explicit
  child species declarations are retained; assays and alignment are unchanged.

# SpaMTPData 0.99.4

* Registered native resource release 1.1.0 from published Zenodo record
  22733262. The seven SpatialExperiment files have their own immutable URLs,
  sizes and MD5 values and are now selected by the default latest version.
* Retained a complete 18-resource snapshot, including eleven unchanged
  resources, plus all historical 1.0.0 entries and their original checksums.
* Updated resource documentation and regression tests to distinguish native
  defaults from explicitly requested historical Seurat archives. Publication
  and package metadata do not imply completed ExperimentHub ingestion.

# SpaMTPData 0.99.3

* Added source attribution and a maintainer registration recipe for native
  files in the Zenodo resource family jointly maintained through SpaMTPdb.
  SpaMTPData remains the experimental-data interface and gains no runtime
  dependency on SpaMTPdb, SpaMTP or Seurat.
* Registration requires a published record, verified checksums and fresh
  independent SpatialExperiment reads. It generates candidate registry/Hub
  metadata, retaining historical entries and a complete 18-resource snapshot.
* Published 1.0.0 entries remain unchanged until a new native release is public
  and its candidate metadata has been reviewed; no draft URLs are advertised.

# SpaMTPData 0.99.2

* Added explicit roxygen2 importFrom declarations for every Imports package,
  including S4 replacement accessors, and regression coverage for namespace
  dependency declarations.
* Added an explicit conversion plan and maintainer recipe for the seven
  historical Seurat archives. The recipe verifies original files and saves
  portable SpatialExperiment RDS with counts/data, feature/pixel annotations,
  source-image centroids, paired transcriptomes and conversion provenance.
* Native preparation omits Seurat graphs, reductions, commands, weights,
  optical rasters and polygon/molecule payloads. It runs once with optional
  Seurat support; using the resulting RDS does not require Seurat.
* Staging version 1.1.0 has its own manifest and checksums. Preparing these
  files does not overwrite the published version 1.0.0 registry or upload data.

# SpaMTPData 0.99.1

* Replaced PascalCase resource exports with spaMTPData(), spaMTPDataResource(),
  spaMTPDataResources() and spaMTPDataVersion(), coordinated with SpaMTP >= 0.99.3.
* Added spaMTPExampleData(): a small deterministic synthetic SpatialExperiment
  with paired transcriptomics, requiring neither Seurat nor downloads.
* Preserve published resource 1.0.0 files, classes and hashes. Registry rows
  identify archived Seurat objects with requires_conversion; no implicit
  conversion or relabelling as native resources is performed.
* Moved the resource manifest to inst/manifest and added consistency checks
  against Hub metadata. Only Hub metadata remains under inst/extdata.
* Added local-file size/checksum and object-class checks, versioned download
  caches, and offline reuse of verified cached files. A development-only
  verify=FALSE bypasses local checksums but never class validation.
* Updated the vignette and documentation for the coordinated native workflow.

# SpaMTPData 0.99.0

* Initial ExperimentHub package scaffold for SpaMTP demonstration data.
* Registered existing Zenodo-hosted mouse brain, bladder, simulated
  single-cell, annotation-refinement, import-example and human brain resources.
* Added local-resource support for offline development and testing.
* Added a verified, persistent Zenodo fallback for the pre-ingestion period.
