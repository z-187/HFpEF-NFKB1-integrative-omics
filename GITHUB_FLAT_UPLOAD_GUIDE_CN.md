# GitHub 网页上传：必须按此方法操作

## 关键原因

不要从 ZIP 压缩包窗口直接拖文件到 GitHub。  
从压缩包窗口拖动会导致文件被打散、自动改名，甚至文件名与内容错位。

## 正确步骤

1. 下载 `HFpEF-NFKB1-integrative-omics_FLAT_WEB_UPLOAD.zip`。
2. 在 Windows 中右键该 ZIP：
   **提取全部（Extract All）**。
3. 进入正常的黄色文件夹：
   `HFpEF-NFKB1-integrative-omics_FLAT_UPLOAD`
4. 确认地址栏中没有 `.zip`，文件夹类型不是“压缩文件夹”。
5. 在该正常文件夹中按 `Ctrl+A` 全选所有文件。
6. 拖到 GitHub 的 Upload files 页面。
7. 上传列表中所有文件都应是唯一名称；不应出现：
   - `README (1).md`
   - `README (2).md`
   - `download`
8. 提交说明填写：
   `Add complete flat-layout reproducible workflow for OMI-2026-0142 revision`

## 上传后必须核对

- GitHub 首页显示 `README.md`
- `CITATION.cff` 第一行是 `cff-version: 1.2.0`
- `LICENSE` 第一行是 `MIT License`
- `06_stage6_TF_dependent_macrophage_vascular_communication_FINAL_v3.R` 开头显示 `Stage 6 FINAL v3`
- `01b_stage1_metadata_patch_FIXED_v1.R` 存在
- 不存在 `README (1).md`、`download` 或任何密码文件
