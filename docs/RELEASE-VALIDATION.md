# Release candidate validation — 2026-10-07

Base commit: `3a610f0b13541194cbd4667e02fc22f5a7f1f512`, plus the release-preparation changes in this checkout.

## Results

- Platform: macOS Apple Silicon; R 4.5.2; OlinkAnalyze 5.0.2; testthat 3.3.2.
- `R CMD build .`: passed, including both HTML vignettes.
- `R CMD check --no-manual --as-cran NPXplore_1.0.0.tar.gz`: 0 errors,
  0 warnings, 2 notes. The incoming check identifies a new submission;
  the timestamp check could not verify the current time. URL checks were
  repeated with network access; no invalid-URL findings remained.
- Installed-package tests: **396 passes, 0 failures, 0 warnings, 0 skips**.
  The RStudio-bundled Pandoc was added to PATH for HTML export tests.
- `Rscript tools/check-examples.R`: passed on the real checkout demo inputs.
  Condition (Cancer minus Healthy): 5,415 assay results, 1,157 DEPs.
  Group (ANOVA): 5,415 assay results, 2,634 DEPs.
- The checked archive installs as 1.0.0, contains the default annotation,
  and exposes both HTML vignettes.
- Document input paths and workflow YAML syntax checked.

Archive SHA-256: `fa24f8ef290fae86cce71bb29b7b06351bfdb51d65837b53122b460bde3d2ffb`.

## Verification limits before publishing

- No GitHub Actions run has been triggered by these local changes. The new
  workflow covers Linux R 4.5.0 and current R on Linux/macOS/Windows; those
  environments are pending, not locally verified.
- The test suite ran with testthat 3.3.2. Its declared 3.2.1 floor corresponds
  to support for mocking non-function objects; that exact old version was not
  executed locally.
- The supported R floor is 4.5.0; local verification used 4.5.2. OlinkAnalyze
  5.0.2 is the checked dependency baseline, not a proven historical minimum.
- Full live Enrichr/STRING analysis and thousands of default plot exports were
  not run. Offline tests exercise mocked network responses and the demo smoke
  check exercises QC/statistics, not service availability.
- PDF manual generation was not checked (`--no-manual`). HTML vignettes were
  built and checked.
- No release tag or public GitHub Release has been created. Follow
  [RELEASING.md](RELEASING.md) after remote checks pass.
