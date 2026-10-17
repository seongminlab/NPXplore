# NPXplore 1.0.0

* Read Olink Explore HT NPX data, join metadata, and report sample/assay QC.
* Run independent or paired two-group tests, ANOVA and post-hoc analysis,
  DEP selection, and multi-group clustering.
* Export tables, static and interactive visualizations, GO/KEGG enrichment,
  and STRING protein interaction networks with stage diagnostics.
* Correct checkout example paths and document optional vignette installation.
* Regenerate reference exports with default QC and current input metadata;
  exclude final FAIL samples and retain WARN samples. Validate metadata,
  NPX sample columns, statistical tables, and DEP membership against fresh
  default analyses (Condition: 1,157 DEPs; Group: 2,634 DEPs).
* Add reproducible export/check scripts with Git base commit, working-tree
  state, source/input SHA-256 hashes, parameters, and session information.

## Example data provenance

The examples use publicly deposited experimental Olink Explore HT NPX measurements from [PRIDE accession PAD000002](https://www.ebi.ac.uk/pride/archive/projects/PAD000002), associated with [Beimers et al., *Technical Evaluation of Plasma Proteomics Technologies*](https://doi.org/10.1021/acs.jproteome.5c00221). Repository-supplied metadata select the healthy/cancer and plasma dilution subsets for demonstrating NPXplore. The exported results reflect NPXplore’s documented inputs, QC rules, and statistical settings; they constitute a separate reanalysis from the original publication.
