# Data and resource manifest

Raw expression matrices and large intermediate objects are not redistributed. Download the public resources under their original repository terms and place them under the separate local project root defined by `HFPEF_PROJECT_DIR`.

| Resource | Role | Boundary |
|---|---|---|
| GSE237156 | Stage 2 pharmacotranscriptomic discovery | Sorted CCR2+ and CCR2− macrophages; two replicates per experimental cell. |
| GSE236585 | Primary cardiac single-cell context | Biological sample, not cell, is the inferential unit. |
| GSE236584 | Matched whole-heart bulk support | Same study as GSE236585; not an independent cohort. |
| GSE208425 | Related metabolic immune-cell context | Contextual support only. |
| GSE245034 | External bulk treatment evaluation | Same study system as GSE249412. |
| GSE249412 | Cell-type-resolved treatment evaluation | Not counted as a second independent cohort relative to GSE245034. |
| SCP3342 | Independent human myocardial single-nucleus evaluation | 19 HFpEF and 24 control donors. |
| NicheNet-v2 mouse resources | Stage 6 ligand–target and ligand–receptor prior knowledge | Fixed Zenodo record 7074291; resource hashes are checked in the script. |

Expected filenames, sample inclusion rules, contrasts, parameters, and random seeds are encoded in the staged scripts and manuscript Supplementary Material.
