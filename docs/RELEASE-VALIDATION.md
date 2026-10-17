# Revision validation — 2026-10-16

Package version: **NPXplore 1.0.0**. Both full default example runs executed the clean analysis-source snapshot **`bd5eb80baec7edfdc9ac7d9bb7f7429add79c87e`**, based on `f7e3f4e1c5a3cd7d2e29fcffd8c22ac002d1efa2`. The snapshot is retained in `source-revision.bundle` in NPXplore-examples. Original working branches were not committed or moved. This validation report was added after the analysis snapshot; it does not change the executed analysis code.

## Inputs and default results

The examples use publicly deposited experimental Olink Explore HT NPX measurements from [PRIDE accession PAD000002](https://www.ebi.ac.uk/pride/archive/projects/PAD000002), associated with [Beimers et al., *Technical Evaluation of Plasma Proteomics Technologies*](https://doi.org/10.1021/acs.jproteome.5c00221). Repository-supplied metadata select the healthy/cancer and plasma dilution subsets for demonstrating NPXplore. The exported results reflect NPXplore’s documented inputs, QC rules, and statistical settings; they constitute a separate reanalysis from the original publication.

The package and examples copies of the parquet and both metadata inputs are byte-identical. Measurements SHA-256: `89d86042905610bb449af112ea5d3947f9b6d27c4b945c41bd11f225c21e33bf`.

| Example | Input metadata rows | Retained samples | Tested assays | DEPs |
| --- | ---: | ---: | ---: | ---: |
| Condition, independent t-test | 43 | 42 (22 Healthy, 20 Cancer) | 5,415 | 1,157 |
| Group, ANOVA | 30 | 30 (3 per group) | 5,415 | 2,634 |

These are fresh complete pipeline runs into empty output directories, using `npx_pipeline()` defaults. Only data, metadata and output_dir were supplied, plus `variable = "Group"` for the multi-group analysis. Defaults include `include_warn = TRUE`, `include_fail = FALSE`, `logFC = 0`, `adj.p = 0.05`, all assay boxplots, raw-NPX clustering, and live functional analysis. The driver records its outer random seed; UMAP, clustering and network functions also use their recorded source defaults.

O5-25 (final FAIL) is absent from analysis metadata and the NPX matrix. O5-10, O5-21 and O5-26 (final WARN) are retained. Original QC and raw-data exports intentionally preserve excluded samples. Input IDs, repeat IDs and group labels are preserved; the multi-group output uses Control and Dilute_1 through Dilute_9. The previous 1,241-DEP export and old concentration labels are superseded.

## Fresh validation evidence

- Full local package tests: **406 passes; 0 failures, errors, warnings or skips**. This includes actual Wilcoxon execution, the public analysis API, and pipeline FAIL filtering with metadata/NPX export consistency.
- `R CMD build`: passed, including both HTML vignettes.
- `R CMD check --no-manual --as-cran`: **0 errors, 0 warnings, 1 NOTE**. The NOTE is `unable to verify current time`. Remote CRAN incoming checks were disabled; this is local package validation.
- Installed package: 1.0.0; public analysis functions, annotation file and both HTML vignettes present.
- `Rscript tools/check-examples.R <regenerated-output-root>`: passed for both complete exports. The check recomputes default QC and analyses and compares both metadata files, every QC CSV, raw NPX, the saved NPX matrix, statistics, DEP membership and ANOVA posthoc values. It also checks exact source/input hashes, defaults, sample membership, counts and recorded stage statuses.
- The checker rejects the previous examples at `metadata_used.csv`; failed-run manifest regression records Failed, the error and a completion timestamp.
- Recovered the source bundle in a separate checkout, restored the supplied input bytes and verified all 14 source/input file hashes against the recorded manifests.
- New UMAP figures inspected: current Condition and Group labels appear as supplied.

Checked source archive SHA-256: `0d3b129e30397328774ec37294a73df6c73a2801cb2d0fa2318d0a351b7235bb`.

## Complete pipeline stage results

| Example | Assay boxplots | Clustering | GO | PPI |
| --- | --- | --- | --- | --- |
| Paired_sample_test | Success | NotRun | Success | Success |
| multiple_group_sample_test | Success | Success | Success | Success |

Each example's `log/run_manifest.csv` records source commit, version, clean-source status, hashes, seed, defaults, UTC start/completion, result counts and stage outcomes. `source_manifest.csv`, `input_manifest.csv`, `pipeline_defaults.csv`, `sample_membership.csv` and service-specific diagnostics supply the detailed evidence. `sessionInfo.txt` records the R and package versions. The driver finished at 2026-10-16 09:39:38 UTC.

## Reproduction and limits

Use the NPXplore-examples README and source bundle to recover this exact analysis source. Copy its original parquet and metadata CSV files into the checkout before checking the stored exports, so CSV line endings match the input hashes. Run:

```sh
Rscript tools/regenerate-examples.R ../NPXplore-reference
Rscript tools/check-examples.R ../NPXplore-reference
```

R and dependency versions: R 4.5.2, OlinkAnalyze 5.0.2, testthat 3.3.2. Package versions are recorded in each example’s session log. Enrichr and STRING responses can change as databases are updated; their observed statuses above apply to this run. The exact minimum dependency versions and the PDF manual were not covered by these checks. No release tag or GitHub publication was performed.
