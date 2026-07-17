# Run order and script roles

The workflow is staged and stateful. Each primary script validates the required upstream outputs. Run scripts from a fresh R session unless explicitly stated otherwise.

## Primary clean-start chain

1. Stage 1 data audit (`01a`)
2. Stage 1 metadata patch (`01b`)
3. Stage 1 metadata lock (`01c`)
4. Stage 2 GSE237156 directionally opposed discovery
5. Stage 3 initial GSE236585 processing (`03a`)
6. Stage 3 final projection/annotation patch (`03b`)
7. Stage 4 macrophage TF/regulon analysis
8. Stage 5 network-constrained multi-TF sensitivity analysis (`05a`)
9. Stage 5B bootstrap and matched-null analysis (`05b`)
10. Stage 6 complete NicheNet communication analysis
11. Stage 7 descriptive sample-level ridge attribution
12. Stage 8 multicohort and human myocardial evaluation
13. Revision benchmarking and ablation
14. Final figures, source data and supplementary tables

## Recovery utilities

- `R/recovery/03c_stage3_v4_disk_resume_legacy.R` resumes a historical downstream delivery section from an existing checkpoint. It is not required for a clean run of `03b`.
- `R/recovery/06b_stage6_delivery_repair_v1.R` repairs Stage 6 delivery/audit outputs from an existing scientific checkpoint. It does not perform NicheNet analysis and must not replace the complete Stage 6 script.
