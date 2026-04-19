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

<img width="495" height="446" alt="droplet plot" src="https://github.com/user-attachments/assets/9f4d810e-6fae-45ea-9b44-a2d4c0b61c06" />
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
# Advanced scRNA-seq Analysis Tutorial



---

## Overview

This tutorial covers three advanced topics in single-cell RNA-seq analysis:

1. **Complete scRNA-seq Workflow** — From raw counts to cell type identification
2. **Gene-Level Analysis** — Functional enrichment and genomic ranges
3. **Gene Set Enrichment Analysis (GSEA)** — Rank-based pathway analysis

---

## Part 1: Complete scRNA-seq Workflow

### 1.1 Setup and Data Loading

**Steps:**
- Import necessary libraries (`scanpy`, `anndata`, `pooch`)
- Set consistent figure style
- Download sample dataset (two 10x Genomics samples: `s1d1` and `s1d3`)
- Load as `AnnData` objects and merge

---

### 1.2 Quality Control (QC)

**Gene identification patterns:**

| Pattern | Meaning |
|---------|---------|
| `MT-` | Mitochondrial genes |
| `RPS`, `RPL` | Ribosomal genes |
| `^HB[^(P)]` | Hemoglobin genes |

**Quality metrics calculated:**
- Number of genes per cell
- Total counts per cell
- Percentage of mitochondrial counts

**Visualization:**
- Violin plots for metric distributions
- Scatter plots (total counts vs. genes, colored by mitochondrial %)

**Filtering:**
- Remove cells with very few genes
- Remove genes found in very few cells

---

### 1.3 Normalization and Feature Selection

| Step | Method | Outcome |
|------|--------|---------|
| Normalization | Median total counts + log-transform | Comparable counts |
| Feature selection | Identify top 2,000 highly variable genes | Genes that distinguish cell types |

---

### 1.4 Dimensionality Reduction and Clustering

| Method | Purpose |
|--------|---------|
| PCA | Linear transformation to reduce complexity |
| UMAP | 2D projection for visualization |
| Leiden | Graph-based clustering |

**Clustering resolutions to explore:**
- 0.02 — fewer, broader clusters
- 0.5 — balanced clustering
- 2 — more, finer clusters

---

### 1.5 Doublet Detection and Removal

**Tool:** `scrublet`

**Process:**
1. Predict doublets based on gene expression profiles
2. Identify cells with high doublet scores
3. Remove predicted doublets
4. Replot UMAP to visualize cleaned data

---

### 1.6 Cell Type Annotation

**Immune cell marker genes:**

| Cell Type | Marker Genes |
|-----------|---------------|
| CD14+ Monocytes | CD14, LYZ |
| NK cells | GNLY, NKG7 |
| CD4+ T cells | CD4, IL7R |
| CD8+ T cells | CD8A, CD8B |
| B cells | CD79A, MS4A1 |

**Methods:**
- Dot plots for marker gene expression across clusters
- Manual annotation based on marker patterns
- `celltypist` for automated annotation

---

## Part 2: Gene-Level Analysis

### 2.1 Introduction to Gene Ontology (GO)

**Gene Ontology categories:**

| Category | Description |
|----------|-------------|
| Biological Process (BP) | Biological objectives (e.g., cell cycle) |
| Molecular Function (MF) | Molecular activities (e.g., kinase activity) |
| Cellular Component (CC) | Locations (e.g., nucleus, membrane) |

**Database used:** `org.Hs.eg.db` (human gene annotations)

---

### 2.2 Differential Expression and Gene List Creation

**Input gene list examples:**
- "SYNGAP1 interactors"
- Genes from a specific study

**Background sets:**
- All genes detected in transcriptome
- All known protein-coding genes (gene universe)

---

### 2.3 Performing Enrichment Analysis

**Tool:** `clusterProfiler`

**Method:** Over-representation test

**Output:** Functional categories (GO terms or KEGG pathways) significantly over-represented in the input gene list.

---

### 2.4 Visualizing Enrichment Results

| Plot Type | Purpose |
|-----------|---------|
| Bar plot | Most significantly enriched terms |
| Dot plot | Enrichment statistics (GeneRatio, p-value) |
| Network plot (enrichment map) | Relationships between terms |
| Category-gene-network | Which genes drive enrichment |
| Heatmap | Comparison across multiple gene lists |

---

### 2.5 Genomic Ranges Analysis

**Tools:** `GenomicRanges` package

**Input data:**
- Mouse recombination hotspots
- PRDM9 binding sites

**Analysis:** Find overlaps between two sets of genomic ranges

**Output:** Identify which PRDM9 binding sites fall within recombination hotspots

**Visualization:** Plot the overlap

---

## Part 3: Gene Set Enrichment Analysis (GSEA)

### 3.1 What is GSEA?

**Key difference from over-representation tests:**

| Method | Input | Approach |
|--------|-------|----------|
| Over-representation | Gene list (cutoff-based) | Counts genes in sets |
| GSEA | Ranked list of ALL genes | Considers gene order |

**Advantage of GSEA:** Doesn't require arbitrary cutoffs for significance.

---

### 3.2 Setup and Data Preparation

**Input data:** Pre-calculated differential expression table with log2 fold-change for all genes

**Ranked list creation:**
1. Sort genes from most upregulated to most downregulated
2. Convert gene symbols to Entrez IDs

---

### 3.3 Performing GSEA

**Tool:** `clusterProfiler::GSEA()`

**Gene set database:** Molecular Signatures Database (MSigDB)

**Output metric:** Normalized Enrichment Score (NES)

| NES Value | Interpretation |
|-----------|----------------|
| Positive > 0 | Genes enriched at top (upregulated) |
| Negative < 0 | Genes enriched at bottom (downregulated) |
| |NES| > 1 | Significant enrichment |

---

### 3.4 Interpreting GSEA Results

**Visualizations:**

| Plot | Shows |
|------|-------|
| Dot plot | Most significant pathways with NES |
| Classic GSEA plot | Running enrichment score with peak at enrichment |
| Heatmap | Comparison across multiple ranked lists |

**Example pathway:** "Epithelial Mesenchymal Transition" (EMT)

**What a sharp peak indicates:** Pathway genes are highly enriched at the top of the ranked list (upregulated).

---

## Summary of Methods

| Analysis | Purpose | Key Tools |
|----------|---------|-----------|
| QC & Filtering | Remove low-quality cells/genes | scanpy, scrublet |
| Normalization | Make counts comparable | Log-normalization, CPM |
| Dimensionality reduction | Visualize high-dim data | PCA, UMAP, t-SNE |
| Clustering | Identify cell groups | Leiden |
| Annotation | Assign cell types | Marker genes, celltypist |
| Enrichment | Find functional themes | clusterProfiler |
| Genomic ranges | Find spatial overlaps | GenomicRanges |
| GSEA | Rank-based pathway analysis | clusterProfiler::GSEA |

---

## Next Steps

After completing this tutorial, you can:
- Apply these methods to your own scRNA-seq data
- Integrate multiple datasets
- Perform trajectory inference
- Identify rare cell populations
- Publish reproducible analyses

---
# AnnData & ScanPy Advanced Tutorial

**Author:** Aqba Ejaz

---

## Overview

This tutorial provides an in-depth exploration of **AnnData**, the core data structure for single-cell genomics, and **ScanPy** for single-cell data processing. Building on introductory material, this guide focuses on practical data manipulation, normalization, visualization, and understanding relationships between different AnnData components.

---

## Table of Contents

1. [Environment Setup](#1-environment-setup)
2. [Loading a Real Dataset](#2-loading-a-real-dataset)
3. [Understanding the AnnData Structure](#3-understanding-the-anndata-structure)
4. [Exploring Sparse Matrix Storage](#4-exploring-sparse-matrix-storage)
5. [Working with Layers](#5-working-with-layers)
6. [Visualizing Normalization Effects](#6-visualizing-normalization-effects)
7. [Exploring Cell Metadata (.obs)](#7-exploring-cell-metadata-obs)
8. [Adding New Metadata Columns](#8-adding-new-metadata-columns)
9. [Exploring Gene Metadata (.var)](#9-exploring-gene-metadata-var)
10. [Renaming Observations and Variables](#10-renaming-observation-and-variable-names)
11. [Subsetting Data](#11-subsetting-data)
12. [Filtering Cells Based on Quality](#12-filtering-cells-based-on-quality)
13. [Exploring Unstructured Data (.uns)](#13-exploring-unstructured-data-uns)
14. [Working with Multi-dimensional Data (.obsm)](#14-working-with-multi-dimensional-data-obsm)
15. [Exploring Pairwise Data (.obsp)](#15-exploring-pairwise-data-obsp)
16. [Understanding Views vs. Copies](#16-understanding-views-vs-copies)

---

## 1. Environment Setup

### Required Libraries
- `anndata` — core data structure
- `scanpy` — single-cell analysis
- `matplotlib` — visualization
- `numpy` — numerical operations
- `pooch` — data download

### What happens
The system downloads and installs all dependencies for handling single-cell data, including array-api-compat and zarr.

---

## 2. Loading a Real Dataset

### Dataset Information
| Property | Value |
|----------|-------|
| Dataset | PBMC 3k (pre-processed) |
| Cells (observations) | 2,638 |
| Genes (variables) | 11,505 |
| Cell type | Peripheral blood mononuclear cells (PBMCs) |

### Steps
1. Use `pooch` library to download the dataset
2. Read the H5AD file into an AnnData object

---

## 3. Understanding the AnnData Structure

| Slot | Purpose | Content |
|------|---------|---------|
| `.obs` | Cell metadata | n_genes, percent_mito, n_counts, cell type labels |
| `.var` | Gene metadata | gene_names, n_cells, gene_ids |
| `.uns` | Unstructured data | Clustering parameters, colors, PCA results |
| `.obsm` | Multi-dimensional embeddings | PCA, t-SNE, UMAP coordinates |
| `.layers` | Alternative matrices | Raw counts |
| `.obsp` | Pairwise matrices | Distance matrix |

---

## 4. Exploring Sparse Matrix Storage

### Key Insight
Single-cell data is very sparse. In this dataset, only **~6.8%** of entries are non-zero.

### What to examine
- The `.X` matrix structure
- Non-zero values and their indices
- Fraction of non-zero entries

**Why it matters:** Sparse storage is essential for memory efficiency.

---

## 5. Working with Layers

### Steps
1. Copy raw counts from "raw" layer to new layer "counts_per_million"
2. Normalize to counts-per-million (CPM) using `sc.pp.normalize_total`
3. Compare raw and normalized counts for specific genes

### Why use layers
Layers allow storing multiple versions of data without overwriting the original `.X` matrix.

---

## 6. Visualizing Normalization Effects

### Visualization
Create a side-by-side matrix plot comparing:

| Plot | Data |
|------|------|
| Left | CPM-normalized data |
| Right | Raw counts |

### Marker Genes Used
| Gene | Cell Type Marker |
|------|------------------|
| CD8A | CD8+ T cells |
| CD4 | CD4+ T cells |
| KLRB1 | NK cells |

### What you'll observe
CPM normalization adjusts for sequencing depth differences, making expression levels comparable across cells.

---

## 7. Exploring Cell Metadata (.obs)

### Key Columns

| Column | Description |
|--------|-------------|
| n_genes | Number of genes detected per cell |
| percent_mito | Percentage of mitochondrial reads |
| n_counts | Total UMI counts per cell |
| louvain_cell_types | Cluster assignments |

### Tasks
- View the entire `.obs` DataFrame
- List all column names
- Examine cell type annotations as categorical variable
- Count B cells (should be ~342)

---

## 8. Adding New Metadata Columns

### Quality Control Filter
Create a boolean column `is_low_quality` based on:
is_low_quality = percent_mito > 0.03 (3%)

text

**Best practice:** High mitochondrial content (>3-5%) often indicates stressed or dying cells.

---

## 9. Exploring Gene Metadata (.var)

### Key Columns

| Column | Description |
|--------|-------------|
| gene_names | Gene symbols (e.g., CD8A, CD4) |
| n_cells | Number of cells expressing the gene |
| gene_ids | ENSEMBL identifiers (e.g., ENSG000001...) |

---

## 10. Renaming Observation and Variable Names

### Why this matters
Different analyses require different identifier types:

| Identifier Type | Use Case |
|----------------|----------|
| Gene symbols | Readability, visualization |
| ENSEMBL IDs | Stability, database queries |

### Tasks
1. View current observation names (cell barcodes)
2. Change variable names to ENSEMBL IDs
3. Change back to gene symbols

---

## 11. Subsetting Data

### Example Subset
- **Cells:** First 5 cells
- **Genes:** LYZ, FOS, MALAT1
- **Result shape:** 5 × 3

### What to examine
- Log-normalized counts (`.X`)
- Raw counts (`.layers["raw"]`)
- Cell annotations for the subset

---

## 12. Filtering Cells Based on Quality

### Process
Create new AnnData object containing only high-quality cells:

python
adata_high_quality = adata[~adata.obs["is_low_quality"], :]


## Quality Control Results

After filtering cells based on mitochondrial percentage (>3%), the following changes were observed:

| Dataset | Number of Cells |
|---------|-----------------|
| Original | 2,638 |
| High-quality | 2,257 |

**Removed:** 381 low-quality cells (14.4%)
