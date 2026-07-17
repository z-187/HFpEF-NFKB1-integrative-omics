# HFpEF integrative transcriptomics workflow

Reproducible R code associated with manuscript **OMI-2026-0142**:

> *Integrative Transcriptomics Links Dapagliflozin-Opposed Macrophage Programs to BHLHE40-Associated Regulation and NFKB1-Associated Signaling in HFpEF*

## Repository format

This is a deliberately **flat-layout repository** so that every file can be uploaded reliably through the GitHub web interface. All scripts and documentation are stored in the repository root with unique filenames.

## Scientific scope

The analysis distinguishes:

- **BHLHE40-associated program-level regulatory robustness**
- **NFKB1-associated communication breadth**

It does not establish dapagliflozin-specific causal reversal, direct ligand–receptor signaling, diagnostic performance, or universal cross-model conservation. The discovery-ranked TNF–TNFRSF1A endothelial branch was not directionally conserved in independent human myocardium, whereas PDGFB-related fibroblast branches showed broader external support.

## Primary data resources

- **GSE237156**: pharmacotranscriptomic discovery
- **GSE236585**: primary cardiac single-cell context
- **GSE236584**: matched bulk support from the same study as GSE236585
- **GSE208425**: related metabolic immune-cell context
- **GSE245034 and GSE249412**: two modalities from one external treatment study
- **SCP3342**: independent human myocardial single-nucleus evaluation (19 HFpEF and 24 control donors)

Raw public datasets and large intermediate objects are not redistributed.

## Recommended execution order

1. `01a_stage1_data_audit_FIXED_v2.R`
2. `01b_stage1_metadata_patch_FIXED_v1.R`
3. `01c_stage1_metadata_lock_FIXED_v3.R`
4. `02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2.R`
5. `03a_stage3_GSE236585_scRNA_projection_FIXED_v2.R`
6. `03b_stage3_GSE236585_scRNA_projection_FIXED_v4_PATCH.R`
7. `04_stage4_GSE236585_macrophage_TF_regulon_FIXED_v1.R`
8. `05a_stage5_multiTF_virtual_perturbation_FIXED_v2.R`
9. `05b_stage5B_offline_bootstrap_null_FIXED_v1.R`
10. `06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R`
11. `07_stage7_cross_stage_sample_ridge_attribution_FINAL_v2.R`
12. `08_stage8_multicohort_validation_FINAL_v6.R`
13. `09_revision_benchmark_ablation_FINAL_v3.R`
14. `10_revision_final_figures_tables_FINAL_v3.R`

Recovery-only utilities:

- `03c_stage3_v4_disk_resume_legacy.R`
- `06b_stage6_delivery_repair_v1.R`

The recovery scripts do not replace the corresponding primary analyses.

## Local project root

Before execution, define the local analysis project root:

```r
Sys.setenv(HFPEF_PROJECT_DIR = "D:/HFpEF_project")
```

Optional Stage 8 paths:

```r
Sys.setenv(
  HFPEF_ASCII_PROJECT_LINK = "D:/HFpEF_STAGE8_ASCII_LINK",
  HFPEF_TEMP_DIR = "D:/HFpEF_STAGE8_TEMP"
)
```

## Software environment

The workflow was developed in **R 4.6.0 on Windows**. See:

- `requirements_R.txt`
- `requirements_python.txt`
- `capture_session_info.R`
- `ENVIRONMENT_README.md`

Stage 6 uses fixed official NicheNet-v2 mouse resources from Zenodo record `10.5281/zenodo.7074291` and verifies the published MD5 hashes.

## Reproducibility boundaries

- Biological samples or donors—not individual cells—are inference units.
- Discovery signatures, TF candidates, and communication axes are frozen before external evaluation.
- Large raw matrices, intermediate RDS objects, and generated results are not included.
- Exact methods, parameters, score definitions, random seeds, and sensitivity settings are documented in the scripts and manuscript Supplementary Material.

## License

MIT License. Third-party datasets and resources remain subject to their original terms.

## Citation

Use `CITATION.cff` and cite the associated manuscript. Create a fixed GitHub release after upload and online verification.

## Contact

- Quan Zhou: 18707240983@163.com
- Yu Jiang: 3494106@qq.com
