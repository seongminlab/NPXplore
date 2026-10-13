# Publishing NPXplore 1.0.0

1. Confirm a clean checkout and reconcile `DESCRIPTION`, `NEWS.md`, and README
   version/support declarations. Change the NEWS release-candidate wording to
   the actual release date only when publishing.
2. Install all Imports and Suggests (STRINGdb via Bioconductor). Run
   `testthat::test_local()` and `Rscript tools/check-examples.R` from the root.
   Inspect skips and failures, not just the command exit status.
3. Build with `R CMD build .`, then check the resulting archive with
   `R CMD check --as-cran NPXplore_1.0.0.tar.gz`. Require zero errors and
   warnings and review every NOTE. Install that archive in a temporary library;
   verify the version, bundled annotation and both installed HTML vignettes.
4. Run the full README examples online in separate output directories. Inspect
   `log` and GO/PPI diagnostics; mocked tests do not verify live services.
5. Push the reviewed candidate and require all GitHub Actions matrix checks to
   pass for that exact commit. Preserve check logs and dependency versions.
6. Tag the verified commit `v1.0.0` using an annotated tag, push that tag, and
   create the GitHub Release using NEWS as its release notes. Attach the checked
   `NPXplore_1.0.0.tar.gz` archive, not an unchecked archive from a later commit.
7. Install `seongminlab/NPXplore@v1.0.0` into a fresh library and confirm
   `packageVersion("NPXplore")` is `1.0.0`. With `build_vignettes = TRUE`, confirm
   both vignettes open. Never move an already published release tag.

Local verification results belong in `docs/RELEASE-VALIDATION.md`. A local pass
is not evidence that the GitHub matrix or live external services passed.
