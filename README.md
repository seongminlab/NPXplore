<h1 align="center">NPXplore</h1>

<p align="center">
  <img src="docs/images/NPXplore-logo.png" alt="NPXplore logo" width="300">
</p>

<br>

## Introduction

NPXplore is an R package for analyzing Olink Explore HT protein expression data.
It combines quality control, differential expression testing, visualization, clustering, and functional analysis in one workflow.
Use the complete pipeline for a quick analysis or run public functions individually for greater control.

The examples use publicly deposited experimental Olink Explore HT NPX measurements from [PRIDE accession PAD000002](https://www.ebi.ac.uk/pride/archive/projects/PAD000002), associated with [Beimers et al., *Technical Evaluation of Plasma Proteomics Technologies*](https://doi.org/10.1021/acs.jproteome.5c00221). Repository-supplied metadata select the healthy/cancer and plasma dilution subsets for demonstrating NPXplore. The exported results reflect NPXplore’s documented inputs, QC rules, and statistical settings; they constitute a separate reanalysis from the original publication.

Browse example analysis exports in the [NPXplore-examples repository](https://github.com/seongminlab/NPXplore-examples). Download the example files [here](https://github.com/seongminlab/NPXplore-examples/archive/refs/heads/main.zip).

<br>

## 1. Requirements

NPXplore supports Windows, macOS, and Linux. The supported R series is 4.5 and later (`R >= 4.5.0`), with OlinkAnalyze >= 5.0.2. CI checks R 4.5.0 and the current R release. Earlier R versions are outside this release support policy. Dependency versions are recorded in each check.

**Core dependencies:** `OlinkAnalyze`, `dplyr`, `readxl`, `tibble`, `tidyr`.

**Full analysis and visualization:** `enrichR`, `openxlsx`, `ggplot2`, `ggraph`, `htmlwidgets`, `ggtext`, `scales`, `pheatmap`, `plotly`, `igraph`, `STRINGdb`, `umap`, `visNetwork`, `NbClust`.

**Documentation and testing:** `knitr`, `markdown`, `testthat` (>= 3.2.1).

Install the CRAN dependencies in R:

```r
install.packages(c(
  "OlinkAnalyze", "dplyr", "readxl", "tibble", "tidyr",
  "broom", "car", "emmeans",
  "enrichR", "openxlsx", "ggplot2", "ggraph", "htmlwidgets",
  "ggtext", "scales", "pheatmap", "plotly", "igraph",
  "umap", "visNetwork", "NbClust", "knitr", "markdown", "testthat"
))
```

Install STRINGdb through Bioconductor:

```r
install.packages("BiocManager")
BiocManager::install("STRINGdb", ask = FALSE, update = FALSE)
```

BiocManager selects the Bioconductor release compatible with your R installation. See the [Bioconductor installation guide](https://bioconductor.org/install/).
GO/KEGG analysis requires an internet connection to Enrichr and submits selected protein identifiers. PPI analysis downloads STRING reference data. Self-contained interactive HTML exports may require Pandoc; when unavailable, some exports use accompanying dependency folders.

<br>

## 2. Installation

<br>

The commands below install the development version from `main`. Install the full-analysis dependencies in section 1 before running the complete pipeline, including STRINGdb via Bioconductor.

### Install with pak

```r
install.packages("pak")
pak::pak("seongminlab/NPXplore")
```

<br>

### Install with devtools

Use this as an alternative to pak:

```r
install.packages("devtools")
devtools::install_github("seongminlab/NPXplore")
```


<br>

### Download with git clone

Run in a terminal:

```bash
git clone https://github.com/seongminlab/NPXplore.git
cd NPXplore
```

To install this local checkout, start R from the repository root:

```r
devtools::install_local(".", build_vignettes = TRUE)
```

The `example/` directory is included in the GitHub checkout. It is excluded from the installed R package to avoid bundling the large input file; clone or download the repository to follow the examples.

<br>

## 3. Quick start

Start a fresh R session after installation. Run the following code from the repository root, where the `example/` directory is located.



<br>

### Two groups: Condition

```r
library(NPXplore)

data <- OlinkAnalyze::read_npx("example/NPXfile.parquet")
meta_condition <- read.csv("example/metadata_HC_vs_Cancer.csv")
```
<br>

```r
result_condition <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta_condition,
  variable = "Condition",
  output_dir = "output_Condition"
)
```

`Condition` contains **Healthy** and **Cancer**. The default analysis is an independent t-test, with Log2FC expressed as Cancer minus Healthy.

<br>

### Multiple groups: Group


```r
library(NPXplore)

data <- OlinkAnalyze::read_npx("example/NPXfile.parquet")
meta_group <- read.csv("example/metadata_human_dilute.csv")
```
<br>

```r
result_group <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta_group,
  variable = "Group",
  output_dir = "output_Group"
)
```

`Group` contains **Control**, **Dilute_1**, **Dilute_2**, **Dilute_3**, **Dilute_4**, **Dilute_5**, **Dilute_6**, **Dilute_7**, **Dilute_8**, and **Dilute_9**. The default analysis is ANOVA, followed by post-hoc testing and DEP clustering.

Both examples use the package defaults: retain WARN samples, exclude final FAIL samples, select DEPs with adjusted p-value < 0.05, and run GO/KEGG and STRING PPI. All assay boxplots are saved, which can take considerable time. Console messages report stage starts, outcomes, and final completion.

<br>

### More stringent DEP cutoffs

```r
result_strict <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta_condition,
  variable = "Condition",
  logFC = 1,
  adj.p = 0.01,
  output_dir = "output_Strict"
)
```

The current API applies **absolute Log2FC >= 1** and **adjusted p-value < 0.01**, retaining both increased and decreased proteins. ANOVA has no single directional Log2FC and uses only the adjusted p-value cutoff. GO plot terms retain their separate adjusted p-value <= 0.05 criterion.

<br>

### Wilcoxon test

```r
result_wilcox <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta_condition,
  variable = "Condition",
  test = "wilcox",
  output_dir = "output_Wilcox"
)
```

This runs an independent-sample Wilcoxon test on Healthy and Cancer. For matched samples, provide `pair_id` with a genuine subject identifier present in the metadata; the example metadata does not supply a pairing identifier.

Each output directory contains `metadata_used.csv`, `QCreport_metadata.csv`, `Results/`, `QC/`, `Tables/`, and `log/`. GO outputs are under `Results/GeneOntology`; PPI outputs are under `Results/STRING_PPI`.

<br>

## 4. Analysis pipeline

- [Paired sample test](inst/tutorials/Paired_sample_test.md): a two-group workflow using `Condition`, with independent t-test/Wilcoxon examples and paired alternatives for data with a subject identifier.
- [multiple group sample test](inst/tutorials/multiple_group_sample_test.md): an ANOVA workflow using `Group`, including covariates, matching post-hoc tests, k-means clustering, and cluster-specific functional analysis.

Both tutorials follow the analysis and export sequence: input, QC, NPX tables, statistics, clustering where applicable, UMAP, heatmaps, individual boxplots, GO/KEGG, and PPI.

To install the HTML vignettes as well, install `knitr` and `markdown`, then use:

```r
devtools::install_github("seongminlab/NPXplore", build_vignettes = TRUE)
```

The default GitHub installation does not build vignettes. Open them after the documentation-enabled installation:

```r
vignette("Paired_sample_test", package = "NPXplore")
vignette("multiple_group_sample_test", package = "NPXplore")
```

<br>

## 5. Release validation

From the repository root, after installing the dependencies:

```r
testthat::test_local()
rcmdcheck::rcmdcheck(args = "--as-cran", error_on = "warning")
```

Install the developer tools first with `install.packages(c("rcmdcheck", "pkgload"))`. CI is configured to check the package on Windows, macOS, and Linux and retain check logs. Vignette code is displayed with `eval=FALSE`; building the documents does not execute the analysis. Run the checkout example smoke check separately:

```sh
Rscript tools/check-examples.R
```

This validates both example inputs, QC, and statistical analyses without querying Enrichr or STRING. With the supplied inputs and defaults, the two-group comparison yields **1,157 DEPs** and the multi-group comparison **2,634 DEPs** from 5,415 tested assays each. FAIL samples, including O5-25, are excluded; WARN samples are retained.

Regenerate all reference exports from the current checkout into a new directory, then compare their metadata and numerical results against a fresh default analysis:

```sh
Rscript tools/regenerate-examples.R ../NPXplore-reference
Rscript tools/check-examples.R ../NPXplore-reference
```

The full run saves all assay boxplots and queries Enrichr and STRING. It can take hours. Each run records the package version, Git base commit and working-tree state, source/input SHA-256 hashes, actual default parameters, sample membership, result counts, and stage diagnostics under `log/`. These identify the exact source and inputs even when the checkout has uncommitted revisions. External service results depend on the database versions available at run time.

See [NEWS.md](NEWS.md) for changes and [release validation](docs/RELEASE-VALIDATION.md) for the current inputs, default-result counts, and verification evidence. After `v1.0.0` is published, install that fixed version with:

```r
devtools::install_github("seongminlab/NPXplore@v1.0.0", build_vignettes = TRUE)
```

## License

**AGPL-3.0**. See [LICENSE](LICENSE) for the GNU Affero General Public License, version 3.

<br>

## Cite
