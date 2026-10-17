# Example data

The examples use publicly deposited experimental Olink Explore HT NPX measurements from [PRIDE accession PAD000002](https://www.ebi.ac.uk/pride/archive/projects/PAD000002), associated with [Beimers et al., *Technical Evaluation of Plasma Proteomics Technologies*](https://doi.org/10.1021/acs.jproteome.5c00221). Repository-supplied metadata select the healthy/cancer and plasma dilution subsets for demonstrating NPXplore. The exported results reflect NPXplore’s documented inputs, QC rules, and statistical settings; they constitute a separate reanalysis from the original publication.

Run the README and tutorials from the NPXplore repository root.

- `NPXfile.parquet`: NPX measurements exported by NPX Map; the identical file is named `PAD000002npxfile.parquet` in NPXplore-examples. SHA-256: `89d86042905610bb449af112ea5d3947f9b6d27c4b945c41bd11f225c21e33bf`.
- `metadata_HC_vs_Cancer.csv`: 43 input records (22 Healthy, 21 Cancer), including three repeat-measurement IDs with an `_r` suffix. IDs O5-1 through O5-20 select Healthy samples; O5-21 through O5-40 select Cancer samples. Repeat measurements retain their source SampleID; no pairing identifier is supplied.
- `metadata_human_dilute.csv`: 30 input records, with three samples each in Control and Dilute_1 through Dilute_9. These are plasma dilution samples. Control is the analysis reference label; it does not denote a clinical control cohort. The current generic Dilute labels must be used instead of the reversed concentration labels in older exports.

The source authors' [two-group code](https://github.com/coongroup/Plasma-Tech-Comparison/blob/main/Figures/Figure4.R) and [dilution code](https://github.com/coongroup/Plasma-Tech-Comparison/blob/main/Figures/Figure3.R) describe the original experiment. NPXplore preserves the checked-in metadata rather than claiming to reproduce the original paper's sample selection or statistical procedure.

Input counts precede cleaning and sample QC. The default pipeline retains WARN samples and excludes final FAIL samples, including O5-25. See `log/sample_membership.csv` in regenerated examples for each input sample's final inclusion and QC status. The parquet and metadata inputs are unchanged by regeneration.
