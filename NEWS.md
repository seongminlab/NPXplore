# NPXplore 1.0.0

Initial public release candidate.

* Read Olink Explore HT NPX data, join metadata, and report sample/assay QC.
* Run independent or paired two-group tests, ANOVA and post-hoc analysis,
  DEP selection, and multi-group clustering.
* Export tables, static and interactive visualizations, GO/KEGG enrichment,
  and STRING protein interaction networks with stage diagnostics.
* Correct checkout example paths and document optional vignette installation.
* Support R 4.5 and later (locally checked on 4.5.2) and OlinkAnalyze
  5.0.2 and later; require
  testthat 3.2.1 for tests that mock non-function objects.
* Add cross-platform package checks and an offline example smoke check.

## Limitations

* Demo inputs are artificial and are not suitable for research.
* Large demo inputs are distributed in the GitHub checkout, not the R package.
* GO/KEGG requires Enrichr access and submits selected protein identifiers;
  STRING PPI downloads external reference data. Offline tests use mocks and
  do not establish availability of these services.
* Covariates are supported for ANOVA only; paired tests require a real pairing
  identifier, which is absent from the demonstration metadata.
* Vignette chunks are not evaluated during document builds. Full-pipeline
  exports may be slow; some HTML exports need accompanying dependency folders
  when a self-contained export cannot be produced.
