# Public release checklist

Before publishing or creating a fixed release:

- [ ] Confirm `R/06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R` is present and is the complete NicheNet analysis, not the delivery-repair utility.
- [ ] Confirm `R/06_STAGE6_FULL_SCRIPT_REQUIRED.md` is absent.
- [ ] Confirm no file contains local usernames, passwords, API keys, tokens or credential material.
- [ ] Confirm public scripts use `HFPEF_PROJECT_DIR` rather than a personal absolute path.
- [ ] Parse every R script under the intended R version; full rerunning of large public datasets is separate from syntax validation.
- [ ] Run `environment/capture_session_info.R` in the final local environment when feasible.
- [ ] Verify the author list and repository URL in `CITATION.cff`.
- [ ] Confirm README dataset roles and claim boundaries match the revised manuscript.
- [ ] Confirm `MANIFEST.tsv` hashes match the uploaded files.
- [ ] Verify README, LICENSE, CITATION.cff, Stage 2 and Stage 6 from the online repository after upload.
- [ ] Create GitHub release `v1.0.0` only after the uploaded contents are frozen and verified.
