# Paired sample test

Use NPXplore 1.0.0 or later. Run the main examples from top to bottom. Alternative statistical models are marked explicitly; choose one model before continuing.

This guide compares two groups using `Condition`. The default is an independent-sample test. For matched samples, use the paired examples and supply a real subject identifier.

```r
library(NPXplore)
# Run from the root of your cloned NPXplore repository.
```

## 1. Read NPX file

Read NPX measurements and metadata. The package supplies the default HT annotation.

```r
inputs <- NPXplore::npx_read_files(
  data_file = "example/NPXfile.parquet",
  meta_file = "example/metadata.csv"
)
```

## 2. QC experiment

Calculate QC from the original measurements and save reports under `QC`.

```r
qc <- NPXplore::npx_qc_reports(
  inputs$data, meta = inputs$meta, output_dir = "output_Condition"
)
table(qc$sample_qc$QC)
```

Apply QC and cleaning. This example retains WARN samples and excludes final FAIL samples. The function also saves `metadata_used.csv` and `QCreport_metadata.csv`.

```r
prepared <- NPXplore::npx_apply_qc(
  inputs, qc = qc, include_warn = TRUE, include_fail = FALSE,
  output_dir = "output_Condition"
)
```

Save a QC-colored UMAP using measurements before sample exclusion.

```r
qc_umap <- NPXplore::npx_umap(prepared$qc_data, variable = "QC")
NPXplore::npx_save_plot(qc_umap, "output_Condition/QC/UMAP_QC_status.pdf")
```

Inspect the distributions of the cleaned samples. These functions save their plots automatically.

```r
NPXplore::npx_iqrplot(
  prepared$data, variable = "Condition", output_dir = "output_Condition/Results/Distribution"
)
NPXplore::npx_distribution_boxplot(
  prepared$data, variable = "Condition", output_dir = "output_Condition/Results/Distribution"
)
```

## 3. Save NPX table

Save the original metadata-matched assay data and the cleaned NPX matrix under `Tables`. Only the exported matrix uses assay-median imputation; statistical analysis still uses `prepared$data`.

```r
NPXplore::npx_export_tables(
  prepared$data, raw_data = prepared$raw_data, output_dir = "output_Condition"
)
```

## 4. Statistical analysis

The default DEP cutoff is adjusted p-value < 0.05. Statistical functions return all tested assays in `analysis$result` and selected DEPs in `analysis$deps`.
### Independent t-test

```r
analysis <- NPXplore::npx_ttest(prepared$data, variable = "Condition")
```

### Wilcoxon test: alternative

Use this call instead of the t-test above for an independent-sample Wilcoxon test.

```r
analysis <- NPXplore::npx_wilcox(prepared$data, variable = "Condition")
```

### Paired tests: alternatives for matched samples

Use **one** of the following calls only when metadata contains a subject identifier linking the two measurements. Replace `SubjectID` with that column name; do not use the unique sample ID as the pairing identifier.

```r
analysis <- NPXplore::npx_ttest(
  prepared$data, variable = "Condition", pair_id = "SubjectID"
)

# Alternative: paired Wilcoxon test
analysis <- NPXplore::npx_wilcox(
  prepared$data, variable = "Condition", pair_id = "SubjectID"
)
```

ANOVA and post-hoc tests belong to the multiple-group guide. Covariates are currently supported only for ANOVA.

### Volcano plot

After choosing a test, plot all statistical results with the default DEP thresholds and save the plot.

```r
NPXplore::npx_volcano(analysis, output_dir = "output_Condition/Results/DEPs/VolcanoPlot")
```

### Save statistical results

Save the complete statistical table and DEP table. `save_npx = FALSE` preserves the NPX files already saved in step 3.

```r
NPXplore::npx_export_tables(
  prepared$data, analysis = analysis,
  output_dir = "output_Condition", save_npx = FALSE
)
```

## 5. Clustering analysis

Not required for this two-group workflow. GO separates positive and negative Log2FC DEPs, and PPI uses all DEPs. Continue to UMAP.

## 6. UMAP

Calculate UMAP from the cleaned analysis data and save it.

```r
umap <- NPXplore::npx_umap(prepared$data, variable = "Condition")
NPXplore::npx_save_plot(umap, "output_Condition/Results/UMAP/UMAP.pdf")
NPXplore::npx_save_plot(umap, "output_Condition/Results/UMAP/UMAP.png")
```

## 7. Heatmap

### All assays

Plot NPX values without row-wise standardization, then save the heatmap.

```r
heatmap_all <- NPXplore::npx_heatmap(prepared$data, variable = "Condition", scale_rows = FALSE)
NPXplore::npx_save_plot(heatmap_all, "output_Condition/Results/Heatmap/Heatmap_NPX.pdf")
```

### DEPs

Run this section only if `analysis$deps` contains assays.

```r
heatmap_deps <- NPXplore::npx_heatmap(
  prepared$data, variable = "Condition", assays = analysis$deps$Assay, scale_rows = FALSE
)
NPXplore::npx_save_plot(heatmap_deps, "output_Condition/Results/DEPs/DEPs_heatmap.pdf")
```

### Pathway-specific DEPs

Use the annotation's `Human_Pathway` assignments. This function creates and saves the pathway heatmaps.

```r
NPXplore::npx_pathway_heatmaps(
  prepared$data, analysis, inputs$anno, variable = "Condition",
  output_dir = "output_Condition/Results/DEPs/Heatmap_Pathway"
)
```

## 8. Boxplot for each Assays

Save individual plots for all assays. This can take considerable time on large panels.

```r
NPXplore::npx_dep_boxplot(
  prepared$data, analysis, meta = prepared$meta, variable = "Condition",
  output_dir = "output_Condition/Results/DEPs/Boxplot"
)
```

To save only DEP assays instead, add `assays = analysis$deps$Assay` to the call above.

## 9. Gene Ontology analysis

Analyze positive and negative Log2FC DEPs separately. GO/KEGG queries send selected protein identifiers to Enrichr. The function saves tables and plots under `Results/GeneOntology`.

Select the top 10 terms by Combined Score among terms with adjusted p-value <= 0.05; display those terms in descending GeneRatio order.

```r
functional <- NPXplore::npx_functional_analysis(
  analysis, data = prepared$data,
  output_dir = "output_Condition", padj_cutoff = 0.05, top_n = 10
)
functional$status
functional$diagnostics
```

The diagnostics distinguish empty inputs, no significant terms, server/query errors, and output errors.

## 10. Protein-Protein Interaction

Create and save STRING networks for all DEPs. Pathway networks use the same annotation as the pathway heatmaps.

```r
ppi <- NPXplore::npx_ppi_network(
  analysis, data = prepared$data, anno = inputs$anno,
  output_dir = "output_Condition"
)
ppi$status
ppi$summary
```

All-DEP and pathway networks color nodes by Log2FC.

Networks and node/edge tables are saved under `Results/STRING_PPI`; pathway outputs are under `Results/STRING_PPI/Pathways`. The summary distinguishes insufficient mapping, absent interactions, query failures, and export failures.
