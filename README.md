# NPXplore 1.0.0

<br>

## Introduction

NPXplore is an R package for analyzing Olink Explore HT protein expression data.
It combines quality control, differential expression testing, visualization, clustering, and functional analysis in one workflow.
Use the complete pipeline for a quick analysis or run public functions individually for greater control.

<br>

## 1. Requirements

The package declares R >= 4.1.0; the release was tested with R 4.5.2. Dependencies may require a newer R version than the package minimum.

**Core dependencies:** `OlinkAnalyze`, `dplyr`, `readxl`, `tibble`, `tidyr`.

**Full analysis and visualization:** `enrichR`, `openxlsx`, `ggplot2`, `ggraph`, `htmlwidgets`, `ggtext`, `scales`, `pheatmap`, `plotly`, `igraph`, `STRINGdb`, `umap`, `visNetwork`, `NbClust`.

**Documentation and testing:** `knitr`, `markdown`, `testthat` (>= 3.0.0).

Install the CRAN dependencies in R:

```r
install.packages(c(
  "OlinkAnalyze", "dplyr", "readxl", "tibble", "tidyr",
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

Replace `YOUR_GITHUB_USERNAME` below with the repository owner's GitHub username.

<br>

### Install with remotes

```r
install.packages("remotes")
remotes::install_github(
  "YOUR_GITHUB_USERNAME/NPXplore",
  dependencies = TRUE, build_vignettes = TRUE, upgrade = "never"
)
```

<br>

### Install with devtools

Use this as an alternative to remotes:

```r
install.packages("devtools")
devtools::install_github(
  "YOUR_GITHUB_USERNAME/NPXplore",
  dependencies = TRUE, build_vignettes = TRUE, upgrade = "never"
)
```

The installation examples include suggested packages and build the step-by-step vignettes. See the [remotes installation reference](https://remotes.r-lib.org/reference/install_github.html).

<br>

### Download with git clone

Run in a terminal:

```bash
git clone https://github.com/YOUR_GITHUB_USERNAME/NPXplore.git
cd NPXplore
```

To install this local checkout, start R from the repository root:

```r
remotes::install_local(".", dependencies = TRUE, build_vignettes = TRUE, upgrade = "never")
```

The `example/` directory is included in the GitHub checkout. It is excluded from the installed R package to avoid bundling the large input file; clone or download the repository to follow the examples.

<br>

## 3. Quick start

Start a fresh R session after installation. Run the following code from the repository root, where the `example/` directory is located.

```r
library(NPXplore)

data <- OlinkAnalyze::read_npx("example/NPXfile.parquet")
meta <- read.csv("example/metadata.csv", check.names = FALSE)
```

<br>

### Two groups: Condition

```r
result_condition <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta,
  variable = "Condition",
  output_dir = "output_Condition"
)
```

`Condition` contains **Healthy** and **Disease**. The default analysis is an independent t-test, with Log2FC expressed as Disease minus Healthy.

<br>

### Multiple groups: Group

```r
result_group <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta,
  variable = "Group",
  output_dir = "output_Group"
)
```

`Group` contains **Healthy**, **Group_1**, and **Group_2**. The default analysis is ANOVA, followed by post-hoc testing and DEP clustering.

Both examples use the package defaults: retain WARN samples, exclude final FAIL samples, select DEPs with adjusted p-value < 0.05, and run GO/KEGG and STRING PPI. All assay boxplots are saved, which can take considerable time. Console messages report stage starts, outcomes, and final completion.

<br>

### More stringent DEP cutoffs

```r
result_strict <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta,
  variable = "Condition",
  logFC = 1,
  adj.p = 0.01,
  output_dir = "output_Condition_strict"
)
```

The current API applies **absolute Log2FC >= 1** and **adjusted p-value < 0.01**, retaining both increased and decreased proteins. ANOVA has no single directional Log2FC and uses only the adjusted p-value cutoff. GO plot terms retain their separate adjusted p-value <= 0.05 criterion.

<br>

### Wilcoxon test

```r
result_wilcox <- NPXplore::npx_pipeline(
  data_file = data,
  meta_file = meta,
  variable = "Condition",
  test = "wilcox",
  output_dir = "output_Condition_wilcox"
)
```

This runs an independent-sample Wilcoxon test on Healthy and Disease. For matched samples, provide `pair_id` with a genuine subject identifier present in the metadata; the example metadata does not supply a pairing identifier.

Each output directory contains `metadata_used.csv`, `QCreport_metadata.csv`, `Results/`, `QC/`, `Tables/`, and `log/`. GO outputs are under `Results/GeneOntology`; PPI outputs are under `Results/STRING_PPI`.

<br>

## 4. Analysis pipeline

- [Paired sample test](inst/tutorials/Paired_sample_test.md): a two-group workflow using `Condition`, with independent t-test/Wilcoxon examples and paired alternatives for data with a subject identifier.
- [multiple group sample test](inst/tutorials/multiple_group_sample_test.md): an ANOVA workflow using `Group`, including covariates, matching post-hoc tests, k-means clustering, and cluster-specific functional analysis.

Both tutorials follow the analysis and export sequence: input, QC, NPX tables, statistics, clustering where applicable, UMAP, heatmaps, individual boxplots, GO/KEGG, and PPI.

Open the installed HTML vignettes in R:

```r
vignette("Paired_sample_test", package = "NPXplore")
vignette("multiple_group_sample_test", package = "NPXplore")
```

<br>

## 5. SessionInfo

Inspect your own analysis environment with:

```r
sessionInfo()
```

The release validation environment is recorded in [SessionInfo](docs/SessionInfo.txt). Machine-specific library paths are omitted from that record.

<br>

## License

**AGPL-3.0**. See [LICENSE](LICENSE) for the GNU Affero General Public License, version 3.

<br>

## Cite
