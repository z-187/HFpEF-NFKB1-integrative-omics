# GitHub 网页上传说明

本压缩包已包含完整 Stage 1–10 代码链，包括完整 Stage 2 和 Stage 6 主分析脚本。

## 一、新建仓库

1. GitHub 页面点击 **New**。
2. Repository name：`HFpEF-NFKB1-integrative-omics`
3. Description 建议填写：

   `Reproducible R workflow for integrative transcriptomic analysis of dapagliflozin-opposed macrophage programs, BHLHE40-associated regulation, and NFKB1-associated signaling in HFpEF.`

4. Visibility 选择 **Public**。
5. 不要勾选 Add README、Add .gitignore 或 Add license，因为本包已经包含。
6. 点击 **Create repository**。

## 二、上传

1. 解压 `HFpEF-NFKB1-integrative-omics_FINAL_GITHUB_UPLOAD.zip`。
2. 进入解压后的 `HFpEF-NFKB1-integrative-omics` 文件夹。
3. 在新仓库页面点击 **uploading an existing file**。
4. 将该文件夹中的全部文件和文件夹拖入上传区域。不要上传 ZIP 本身。
5. Commit message：

   `Add complete reproducible workflow for OMI-2026-0142 revision`

6. 点击 **Commit changes**。

## 三、上传后核对

必须确认：

- 首页 README 显示当前 BHLHE40/NFKB1 双层结论；
- `LICENSE` 第一行是 `MIT License`；
- `CITATION.cff` 第一行是 `cff-version: 1.2.0`；
- `R/02_stage2_GSE237156_drug_opposed_discovery_FIXED_v2.R` 存在；
- `R/06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R` 存在；
- `R/recovery/06b_stage6_delivery_repair_v1.R` 存在并明确为 recovery-only；
- 不存在 `06_STAGE6_FULL_SCRIPT_REQUIRED.md`；
- 不存在 `.repository_cleared`、密码文件或 ZIP 文件。

确认无误后，可在右侧 **Releases → Create a new release** 创建 `v1.0.0`。
