# gsd4_trio_wes_2025_minibench

This task reproduces a simplified whole exome sequencing (WES) trio analysis from the paper:

**Clinical phenotype and trio whole exome sequencing data from a patient with glycogen storage disease IV in Indonesia**  
Ivan William Harsono, Yulia Ariani, Beben Benyamin, Fadilah Fadilah, Dwi Ari Pujianto, Cut Nurul Hafifah, Titis Prawitasari  
*Data in Brief* (2025)  
DOI: [10.1016/j.dib.2024.111231](https://doi.org/10.1016/j.dib.2024.111231)

---

## 1. Data sources
- Publicly available sequencing data: **BioProject PRJNA1077459** (SRA accessions: SRR27997290, SRR27997291, SRR27997292).  
- Subsampled to 150,000 read pairs per sample with `seqtk sample -s100`.  
- Reference: **GRCh38 chromosome 3** with `.fai` and `.dict`.  
- Packaged in `GSD4_trio_DataFiles.zip` (≈236 MB).

---

## 2. How to get the data
Data is stored externally and not committed to this repository.  

1. Download `GSD4_trio_DataFiles.zip` from the shared **DataFiles** Drive:  
   ```
   DataFiles/gsd4_trio_wes_2025_minibench/GSD4_trio_DataFiles.zip
   ```
2. Unpack into the repo with:
   ```bash
   bash sample_submission_boilerplate/workflow/00_unpack_data.sh /path/to/GSD4_trio_DataFiles.zip
   ```

This will populate:
```
sample_submission_boilerplate/data/
  fastq/*.mini150k_s100.fastq.gz
  ref/GRCh38_chr3.fa{,.fai,.dict}
```

---

## 3. Workflow steps

All scripts live in `sample_submission_boilerplate/workflow/`.

### Step 0 – Unpack data
Prepare FASTQs and reference locally.  
```bash
bash workflow/00_unpack_data.sh /path/to/GSD4_trio_DataFiles.zip
```

### Step 1 – Alignment
Align reads for father, mother, proband.  
```bash
bash workflow/01_align.sh
```
Outputs: `results/bam/*.mkdup.bam`

### Step 2 – Joint variant calling
Call variants across trio.  
```bash
bash workflow/02_joint_call_bcftools.sh
```
Outputs: `results/vcf/trio.vcf.gz`

### Step 3 – QC and Mendelian check
Compute stats and detect Mendelian violations.  
```bash
bash workflow/03_qc_and_mendel.sh
```
Outputs: `results/qc/variants.bcftools.stats.txt`, `results/qc/mendel.txt`

### Step 4 – Annotation (optional)
Annotate variants and extract candidates in GBE1.  
```bash
bash workflow/04_annotate.sh
```

### Step 5 – Evaluation
Summarize results in JSON.  
```bash
python workflow/05_evaluate.py
```

---

## 4. Questions & Answers

- `sample_submission_boilerplate/questions.yaml`  
- `sample_submission_boilerplate/answers.yaml`

These files define 5 benchmark questions and their ground-truth answers.

---

## 5. Environment

A reproducible conda environment is defined in `workflow/env.yml`:

```bash
mamba env create -f workflow/env.yml
mamba activate gsd4-trio
```

Dependencies include: `bwa-mem2`, `samtools`, `bcftools`, `picard`, `openjdk`, and Python 3.11.

---

## 6. Directory structure

```
sample_submission_boilerplate/
├── answers.yaml
├── data/
│   └── README.md         # explains external DataFiles
├── metadata.yaml
├── questions.yaml
└── workflow/
    ├── 00_unpack_data.sh
    ├── 01_align.sh
    ├── 02_joint_call_bcftools.sh
    ├── 03_qc_and_mendel.sh
    ├── 04_annotate.sh
    ├── 05_evaluate.py
    └── env.yml
```

---

## 7. Runtime expectations
- Data size: ~236 MB (mini FASTQs + chr3 reference).  
- Runtime: ≤ 1 hour on a modest 8-core environment.  
- Human time to reproduce: ~60 minutes.

---
