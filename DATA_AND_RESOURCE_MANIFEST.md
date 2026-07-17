# HFpEF integrative transcriptomics workflow

Reproducible R code associated with manuscript **OMI-2026-0142**:

> *Integrative Transcriptomics Links Dapagliflozin-Opposed Macrophage Programs to BHLHE40-Associated Regulation and NFKB1-Associated Signaling in HFpEF*

## Scientific scope

This repository contains the staged analysis workflow used for the major revision. The analysis separates two complementary candidate layers:

- **BHLHE40**: prioritized mainly by program-level regulatory robustness.
- **NFKB1**: prioritized mainly by the breadth of curated macrophage-to-vascular/stromal communication links.

The workflow is hypothesis-generating. It does **not** establish dapagliflozin-specific causal reversal, a diagnostic classifier, direct ligand-receptor engagement, or universal cross-model conservation. The discovery-ranked TNF-TNFRSF1A endothelial branch was not directionally conserved in independent human myocardium, whereas PDGFB-related fibroblast branches received broader external support.

## Data resources and evidence boundaries

| Resource | Primary role | Boundary |
|---|---|---|
| GSE237156 | Pharmacotranscriptomic discovery in sorted CCR2+ and CCR2− macrophages | The diet and treatment contrasts share the HFD-vehicle group; results are directionally opposed patterns, not an independent estimate of drug-specific reversal. |
| GSE236585 | Primary cardiac single-cell localization and sample-level pseudobulk analysis | Three HFpEF and three control samples; cells are not treated as biological replicates. |
| GSE236584 | Matched whole-heart bulk support | Same study as GSE236585; not an independent cohort. |
| GSE208425 | Related metabolic immune-cell context | Contextual support only. |
| GSE245034 and GSE249412 | External treatment study, bulk and cell-type-resolved modalities | Two modalities from the same study; not two independent cohorts. |
| SCP3342 | Independent human myocardial single-nucleus evaluation | 19 HFpEF and 24 control donors; donor is the unit of inference. |

Raw GEO and Single Cell Portal data are not redistributed. Users must obtain the public datasets under their original repository terms.

## Repository structure

```text
.
├── R/                       # Staged analysis scripts
│   ├── 01a...01c            # Data audit and metadata lock
│   ├── 02                   # GSE237156 directionally opposed discovery
│   ├── 03a...03b            # GSE236585 single-cell projection and final patch
│   ├── 04                   # Macrophage TF/regulon analysis
│   ├── 05a...05b            # Network-constrained sensitivity and bootstrap/null analyses
│   ├── 06                   # Complete NicheNet communication analysis
│   ├── 07                   # Descriptive sample-level ridge attribution
│   ├── 08                   # Multicohort and human myocardial evaluation
│   ├── 09                   # Benchmarking and ablation
│   ├── 10                   # Publication figures, source data and supplementary tables
│   └── recovery/            # Recovery-only utilities; not substitutes for primary scripts
├── config/                  # Environment-variable example
├── data/                    # Data-placement guidance only
├── docs/                    # Run order, data manifest and release guidance
├── environment/             # Software requirements and session-info capture
├── LICENSE
├── CITATION.cff
└── MANIFEST.tsv
```

## Local project layout

The scripts expect a separate local project root containing the downloaded public data and generated stage folders:

```text
<HFPEF_PROJECT_DIR>/
├── 0.GEO/
├── 00_external_resources/
├── 01_stage1_data_audit_FIXED_v2/
├── 01_stage1_metadata_patch_FIXED_v1/
├── 01_stage1_metadata_lock_FIXED_v3/
├── 02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2/
├── 03_stage3_GSE236585_scRNA_projection_FIXED_v2/
├── 03_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH/
├── 04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1/
├── 05_stage5_multiTF_virtual_perturbation_FIXED_v2/
├── 05B_stage5B_OFFLINE_bootstrap_null_FIXED_v1/
├── 06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3/
├── 07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2/
├── 08_stage8_multicohort_validation_FINAL_v6/
└── REVISION_Benchmark_Ablation_FINAL_v3/
```

Set the project root before running the scripts:

```r
Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
```

Optional Windows ASCII-only paths for Stage 8 can be set as follows:

```r
Sys.setenv(
  HFPEF_ASCII_PROJECT_LINK = "D:/HFpEF_STAGE8_ASCII_LINK",
  HFPEF_TEMP_DIR = "D:/HFpEF_STAGE8_TEMP"
)
```

## Recommended execution order

| Order | Script | Role |
|---:|---|---|
| 1 | `R/01a_stage1_data_audit_FIXED_v2.R` | Input-data audit and manifest |
| 2 | `R/01b_stage1_metadata_patch_FIXED_v1.R` | Sample-level metadata repair |
| 3 | `R/01c_stage1_metadata_lock_FIXED_v3.R` | Dataset-specific metadata lock |
| 4 | `R/02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2.R` | Directionally opposed macrophage-program discovery |
| 5 | `R/03a_stage3_GSE236585_scRNA_projection_FIXED_v2.R` | Initial single-cell processing and checkpoint generation |
| 6 | `R/03b_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH.R` | Final doublet correction, annotation, pseudobulk and projection |
| 7 | `R/04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1.R` | TF/regulon inference |
| 8 | `R/05a_stage5_multiTF_virtual_perturbation_FIXED_v2.R` | Network-constrained multi-TF sensitivity analysis |
| 9 | `R/05b_stage5B_offline_bootstrap_null_FIXED_v1.R` | Bootstrap and matched-null robustness analysis |
| 10 | `R/06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R` | Complete NicheNet ligand-receptor-target analysis |
| 11 | `R/07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2.R` | Descriptive sample-level attribution |
| 12 | `R/08_stage8_multicohort_validation_FINAL_v6.R` | Multicohort and human myocardial evaluation |
| 13 | `R/09_revision_benchmark_ablation_FINAL_v3.R` | Benchmarking and ablation |
| 14 | `R/10_revision_final_figures_tables_FINAL_v3.R` | Final figures, source data and supplementary tables |

Recovery utilities in `R/recovery/` operate on existing checkpoints and must not replace the corresponding clean-start scripts.

## Software environment

The analysis was developed under **R 4.6.0 on Windows**. Package requirements are listed in `environment/requirements_R.txt` and `environment/requirements_python.txt`. Run `environment/capture_session_info.R` to record the local installed versions.

Stage 6 uses the fixed official NicheNet-v2 mouse resources from Zenodo record `10.5281/zenodo.7074291` and verifies the resource hashes within the script.

## Reproducibility notes

- Biological samples or human donors, rather than cells, are the units of inference.
- Discovery signatures, TF candidates and communication axes are frozen before external evaluation.
- Large raw matrices, intermediate RDS objects and generated output folders are not included.
- Exact score definitions, thresholds, sensitivity settings, random seeds and software details are documented in the scripts and manuscript supplementary tables.
- `MANIFEST.tsv` records repository file roles, versions, sizes and SHA-256 hashes.

## License

Code is released under the MIT License. Public datasets and third-party resources remain subject to their original licenses and terms.

## Citation

Use `CITATION.cff` and cite the associated manuscript. A fixed GitHub release should be created after the repository contents have been uploaded and verified.

## Contact

- Quan Zhou: 18707240983@163.com
- Yu Jiang: 3494106@qq.com
