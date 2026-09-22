<h1 align="center">NPXplore</h1>

<p align="center">
  <img src="docs/images/NPXplore-logo.png" alt="NPXplore logo" width="300">
</p>

<br>

## Introduction

NPXplore is an R package for analyzing Olink Explore HT protein expression data.
It combines quality control, differential expression testing, visualization, clustering, and functional analysis in one workflow.
Use the complete pipeline for a quick analysis or run public functions individually for greater control.

Browse example analysis exports in the [NPXplore-examples repository](https://github.com/seongminlab/NPXplore-examples). Download the example files [here](https://github.com/seongminlab/NPXplore-examples/archive/refs/heads/main.zip).

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
devtools::install_local(".")
```

The `example/` directory is included in the GitHub checkout. It is excluded from the installed R package to avoid bundling the large input file; clone or download the repository to follow the examples.

> **Caution:** The `example/NPXfile.parquet` file contains artificially generated data and must not be used for research. The information in `example/metadata.csv` is demo data and does not represent real individuals or actual observations.

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

Package: **NPXplore_1.0.0**, R 4.5.2, macOS Tahoe 26.5.2 (Apple Silicon).

[Full environment and package versions](docs/SessionInfo.txt)

<details>
<summary>View environment and package versions</summary>

```text
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.5.2

Matrix products: default
BLAS: Apple Accelerate (library path omitted)
LAPACK: R bundled LAPACK; LAPACK version 3.12.1 (library path omitted)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Asia/Seoul
tzcode source: internal

Base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

NPXplore release:
[1] NPXplore_1.0.0

Dependency versions:
  [1] DBI_1.3.0           bitops_1.1-0        gridExtra_2.3.1     readxl_1.5.0        rlang_1.3.0         magrittr_2.0.5      otel_0.2.0
  [8] compiler_4.5.2      RSQLite_3.53.3      systemfonts_1.3.2   png_0.1-9           vctrs_0.7.3         stringr_1.6.0       pkgconfig_2.0.3
 [15] fastmap_1.2.0       backports_1.5.1     dbplyr_2.6.0        labeling_0.4.3      ggraph_2.2.2        caTools_1.18.4      rmarkdown_2.31
 [22] markdown_2.0        ragg_1.5.2          purrr_1.2.2         bit_4.6.0           xfun_0.60           WriteXLS_6.8.0      cachem_1.1.0
 [29] litedown_0.11       jsonlite_2.0.0      blob_1.3.0          tweenr_2.0.3        broom_1.0.13        R6_2.6.1            stringi_1.8.9
 [36] RColorBrewer_1.1-3  reticulate_1.46.0   car_3.1-5           cellranger_1.1.0    estimability_2.0.0  Rcpp_1.1.2          assertthat_0.2.1
 [43] knitr_1.51          Matrix_1.7-6        igraph_2.3.3        tidyselect_1.2.1    yaml_2.3.12         rstudioapi_0.19.0   abind_1.4-8
 [50] viridis_0.6.5       ggtext_0.1.2        enrichR_3.4         gplots_3.3.0        curl_8.0.0          lattice_0.23-1      tibble_3.3.1
 [57] plyr_1.8.9          withr_3.0.3         S7_0.2.2            askpass_1.2.1       evaluate_1.0.5      polyclip_1.10-7     zip_3.0.2
 [64] xml2_1.6.0          pillar_1.11.1       carData_3.0-6       KernSmooth_2.23-27  plotly_4.12.1       generics_0.1.4      ggplot2_4.0.3
 [71] commonmark_2.0.0    scales_1.4.0        xtable_1.8-8        chron_2.3-63        gtools_3.9.5        glue_1.8.1          pheatmap_1.0.13
 [78] emmeans_2.0.4       tools_4.5.2         data.table_1.18.6.1 RSpectra_0.16-2     openxlsx_4.2.8.1    gsubfn_0.7          visNetwork_2.1.4
 [85] mvtnorm_1.4-2       graphlayouts_1.2.5  tidygraph_1.3.1     grid_4.5.2          plotrix_3.8-14      tidyr_1.3.2         crosstalk_1.2.2
 [92] umap_0.2.10.0       duckdb_1.5.5        ggforce_0.5.0       proto_1.0.0         Formula_1.2-6       cli_3.6.6           textshaping_1.0.5
 [99] NbClust_3.0.1       sqldf_0.4-12        viridisLite_0.4.3   arrow_25.0.1        dplyr_1.2.1         gtable_0.3.6        hash_2.2.6.4
[106] digest_0.6.39       OlinkAnalyze_5.0.2  ggrepel_0.9.8       rjson_0.2.23        STRINGdb_2.22.0     htmlwidgets_1.6.4   farver_2.1.2
[113] memoise_2.0.1       htmltools_0.5.9     lifecycle_1.0.5     httr_1.4.8          gridtext_0.1.6      openssl_2.4.2       bit64_4.8.4
[120] MASS_7.3-66
```

</details>

<br>

## License

**AGPL-3.0**. See [LICENSE](LICENSE) for the GNU Affero General Public License, version 3.

<br>

## Cite
