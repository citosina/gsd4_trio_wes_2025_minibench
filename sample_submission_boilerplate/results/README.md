# Results Directory

This folder is **ignored by git** to prevent large files from being committed.

When you run the workflow, the following subdirectories will be created here:

- `vcf/` → Joint VCF outputs
- `qc/` → QC reports and Mendelian check files
- `annotation/` → Annotated variants (if annotation step is run)

> Reproduce by running the workflow on the subsampled FASTQs and reference in `data/`.
