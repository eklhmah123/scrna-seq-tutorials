# Galaxy scRNA-seq Pipeline: Complete Step-by-Step Guide

**Author:** Aqba Ejaz

---

## Overview

This guide walks you through a complete single-cell RNA-seq analysis using the **Galaxy platform**—a web-based, graphical interface that requires no coding. The pipeline processes 10X Chromium data from raw FASTQ files to a filtered count matrix ready for downstream analysis.

---

## Table of Contents

1. [Data Upload and Organization](#1-data-upload-and-organization)
2. [Understanding 10X Chemistry Versions](#2-understanding-10x-chemistry-versions)
3. [Demultiplexing with RNA STARsolo](#3-demultiplexing-and-quantification-with-rna-starsolo)
4. [Mapping Quality Control](#4-mapping-quality-control)
5. [Producing a High-Quality Count Matrix](#5-producing-a-high-quality-count-matrix)
6. [Final Output](#6-final-output)
7. [Key Takeaways](#key-takeaways)

---

## 1. Data Upload and Organization

### Create a New History
- Open Galaxy and create a new history (project folder)
- Rename it descriptively: `scRNA-seq 10X dataset tutorial`

### Import the Raw FASTQ Data
You will import four FASTQ files:

**Read 1 files (Cell Barcode + UMI):**
- `subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz`
- `subset_pbmc_1k_v3_S1_L002_R1_001.fastq.gz`

**Read 2 files (cDNA insert):**
- `subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz`
- `subset_pbmc_1k_v3_S1_L002_R2_001.fastq.gz`

### Import Annotation Files
- **Gene annotation file:** `Homo_sapiens.GRCh37.75.gtf`
- **Cell barcode whitelist:** `3M-february-2018.txt.gz`

---

## 2. Understanding 10X Chemistry Versions

### Chemistry Background

| Chemistry | Read 2 (cDNA) | Read 1 (CB + UMI) | Total Insert Size |
|-----------|---------------|-------------------|-------------------|
| v2 | 98 bp | 26 bp (16+10) | 124 bp |
| v3 | 91 bp | 28 bp (16+12) | 119 bp |

**Key constants:**
- Cell Barcode (CB) is always **16 bp**
- UMI length increased from 10 bp (v2) to 12 bp (v3)

### Determine Your Data's Chemistry
- Click the galaxy-eye icon on any Read 1 FASTQ file
- Count the basepairs in any read sequence
- **28 bp → v3 chemistry** | **26 bp → v2 chemistry**

---

## 3. Demultiplexing and Quantification with RNA STARsolo

### Purpose
- Matches each read's cell barcode against the whitelist
- Uses UMIs to deduplicate PCR artifacts
- Maps cDNA reads to the reference genome
- Counts UMIs per gene per cell

### RNA STARsolo Parameters

| Setting | Selection |
|---------|-----------|
| Reference genome | Human (Homo Sapiens): hg19 Full |
| Gene model file | Homo_sapiens.GRCh37.75.gtf |
| Type of scRNA-seq | Drop-seq or 10X Chromium |
| Input Type | Separate barcode and cDNA reads |
| Chromium chemistry | v3 |
| UMI deduplication | CellRanger2-4 algorithm |
| Barcode matching | Multiple matches (1MM_multi) |

### Expected Output Files (6 files)
- Log — program information and mapping statistics
- Feature Statistic Summaries — detailed QC metrics
- Alignments (BAM) — mapped reads
- Matrix Gene Counts — count matrix
- Barcodes — detected cell barcode sequences
- Genes — gene names/IDs

---

## 4. Mapping Quality Control

### Using MultiQC for Summary Report
- Run MultiQC tool to aggregate STARsolo log information
- **Key metric:** Uniquely mapped reads percentage (>70% is good)

### Direct Inspection of Feature Statistics

**Key metrics to understand:**

| Category | Metric | Meaning |
|----------|--------|---------|
| Barcode | noNoWLmatch | Reads whose barcode doesn't match whitelist |
| Barcode | yesWLmatchExact | Reads with exact barcode match |
| Barcode | yesCellBarcodes | Number of detected cells |
| Gene | noUnmapped | Reads that didn't map to genome |
| Gene | noNoFeature | Reads mapped but not to any gene |
| Gene | MultiFeature | Reads mapping to multiple genes |

**Expected result:** Approximately **5,200 cells** detected.

---

## 5. Producing a High-Quality Count Matrix

### The Problem
Raw STARsolo output includes many empty droplets → extremely sparse matrix.

### Method 1: DefaultDrops (Cell Ranger emulation)

| Parameter | Setting |
|-----------|---------|
| Input format | Bundled |
| Method | DefaultDrops |
| Expected Number of Cells | 3000 |
| Upper Quantile | 0.99 |
| Lower Proportion | 0.1 |

### Method 2: Barcode Rank Plot (Diagnostic)

**What to observe:**
- **Knee point** — sharp bend in the curve
- **Inflection point** — where curve changes slope

These mark the transition between high-quality cells and empty droplets.

### Method 3: EmptyDrops (Statistical Filtering)

| Parameter | Setting |
|-----------|---------|
| Method | EmptyDrops |
| Lower-bound Threshold | 200 |
| FDR Threshold | 0.01 (1% false positives) |

---

## 6. Final Output

The filtered count matrix contains:

- **Real cells only** (high-quality barcodes)
- **All genes** (still sparse, but manageable)

**Three files in bundled format:**
- `barcodes.tsv` — filtered cell barcodes
- `genes.tsv` — gene names/IDs
- `matrix.mtx` — count matrix
<img width="495" height="446" alt="droplet plot" src="https://github.com/user-attachments/assets/9f4d810e-6fae-45ea-9b44-a2d4c0b61c06" />
<img width="566" height="492" alt="droplet plot2" src="https://github.com/user-attachments/assets/652ffb8f-bf2c-475c-a18a-31fcdafb16ef" />

---

## Key Takeaways

| Concept | Summary |
|---------|---------|
| No coding required | All steps via Galaxy's graphical interface |
| Chemistry matters | v2 vs. v3 determines barcode/UMI location |
| Raw output overestimates cells | Includes many empty droplets |
| Filtering is essential | DropletUtils identifies real cells |
| Multiple methods available | DefaultDrops, Barcode Rank Plot, EmptyDrops |
| Final output ready | Filtered matrix for downstream analysis |

---

## Next Steps

This filtered matrix can now be used with:
- **Scanpy/AnnData** for Python-based analysis
- **Seurat** for R-based analysis
- **RaceID** for rare cell type detection

---

## Author

**Aqba Ejaz**
