# Transcriptomics Pipeline: RNA-Seq Analysis & Differential Expression

Comprehensive pipeline for **RNA-Seq data analysis**, covering read alignment, BAM processing, quality control, and differential expression analysis using **metaseqR2**.

Designed for hands-on training in transcriptomics and suitable for execution on Linux-based systems.

---

## Table of Contents

1. [Overview](#overview)
2. [Workflow Summary](#workflow-summary)
3. [Requirements](#requirements)
4. [Pipeline Steps](#pipeline-steps)

   * [1. Setup](#1-setup)
   * [2. Read Alignment (HISAT2)](#2-read-alignment-hisat2)
   * [3. BAM Processing (SAMtools)](#3-bam-processing-samtools)
   * [4. Differential Expression (metaseqR2)](#4-differential-expression-metaseqr2)
   * [5. Alignment Metrics](#5-alignment-metrics)
5. [Outputs](#outputs)
6. [Quality Control & Interpretation](#quality-control--interpretation)
7. [Troubleshooting](#troubleshooting)
8. [Acknowledgements](#acknowledgements)

---

## Overview

This pipeline processes RNA-Seq data from raw reads to biological insight. It includes:

* Alignment to the human reference genome (**hg19**)
* Conversion and processing of alignment files (SAM → BAM)
* Quality control and statistical analysis
* Differential expression analysis using **metaseqR2**

### Biological Context

The dataset consists of prostatic biopsy transcriptomes from patients with advanced hormone-naive prostate cancer.

**Comparison:**
**Post-docetaxel treatment vs Pre-docetaxel baseline**

Reference: [Tzelepi et al.](https://pubmed.ncbi.nlm.nih.gov/)

---

## Workflow Summary

```bash
FASTQ → HISAT2 → SAM → BAM → Sorting → Indexing → metaseqR2 → DE results
```

---

## Requirements

* Linux/Unix environment
* R ≥ 3.6
* ~10 GB disk space

### Tools

* HISAT2
* SAMtools
* metaseqR2 (Bioconductor)

---

## Pipeline Steps

---

### 1. Setup

```bash
mkdir align && cd align

# Download RNA-Seq data
wget http://epigenomics.fleming.gr/~panos/appbio/human.fastq.gz

# Download reference genome (hg19)
wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
tar -xzf hg19_genome.tar.gz

# Download HISAT2
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2.zip
unzip hisat2.zip

# Download SAMtools
wget https://github.com/samtools/samtools/releases/download/1.23/samtools-1.23.tar.bz2
tar -xf samtools-1.23.tar.bz2
cd samtools-1.23 && ./configure && make && cd ..
```

---

### 2. Read Alignment (HISAT2)

HISAT2 performs **splice-aware alignment**, allowing reads to span exon–exon junctions.

```bash
./hisat2-2.2.1/hisat2 \
  -x hg19_genome/genome \
  -U human.fastq.gz \
  -S human_aligned.sam
```

**Key output:**

* `human_aligned.sam`
* Alignment statistics printed in terminal

---

### 3. BAM Processing (SAMtools)

Convert, sort, and index alignments:

```bash
# Convert SAM → BAM
./samtools-1.23/samtools view -b human_aligned.sam > human_aligned.bam

# Sort BAM
./samtools-1.23/samtools sort human_aligned.bam \
  -o human_aligned_sorted.bam

# Index BAM
./samtools-1.23/samtools index human_aligned_sorted.bam
```

**Outputs:**

* `human_aligned_sorted.bam`
* `human_aligned_sorted.bam.bai`

---

### 4. Differential Expression (metaseqR2)

#### Prepare metadata (`targets.txt`)

```text
samplename    filename              condition    paired    strandedness
pre_1         pre_sample1.bam       pre          TRUE      forward
post_1        post_sample1.bam      post         TRUE      forward
```

---

#### Run analysis in R

```r
library(metaseqR2)

metaseqr2(
  sampleList    = readTargets("targets.txt"),
  contrast      = c("post_vs_pre"),
  normalization = "deseq2",
  statistics    = "deseq2",
  qcPlots       = c(
    "mds","biodetection","countsbio","saturation",
    "correl","boxplot","meandiff","meanvar",
    "volcano","mastat"
  ),
  exportWhere   = "./metaseqr2_output",
  org           = "hg19",
  refdb         = "ensembl",
  report        = TRUE
)
```

**Output:**

* Interactive HTML report
* Differential expression tables
* Normalized counts

---

### 5. Alignment Metrics

```r
library(metaseqR2)

annotation <- loadAnnotation(
  org   = "hg19",
  refdb = "ensembl",
  type  = "gene"
)

metaSeqAlignStats(
  targets   = readTargets("targets.txt"),
  bamRanges = annotation
)
```

---

## Outputs

* `human_aligned_sorted.bam` — processed alignments
* `human_aligned_sorted.bam.bai` — index
* `metaseqr2_output/` — results directory

  * `metaseqr2_report.html`
  * differential expression tables
  * normalized counts

---

## Quality Control & Interpretation

### Key metrics to examine:

* **Alignment rate:** expected >80%
* **MDS plot:** clear separation between conditions
* **Volcano plot:** significantly differentially expressed genes

### Biological interpretation:

* Positive log₂ fold-change → upregulated after treatment
* Negative log₂ fold-change → downregulated after treatment

---

## Troubleshooting

**Low alignment rate (<70%)**

* Check FASTQ quality
* Verify reference genome

**metaseqR2 errors**

* Ensure correct paths in `targets.txt`
* Check working directory in R (`getwd()`)

**Memory issues**

* Ensure sufficient RAM (>8 GB recommended)


---

## Acknowledgements

- MSc Applied Bioinformatics (AUTH / IHU)  
- Course: *Proteomics and Functional Genomics*  
- The metaseqR2 pipeline was developed by Dr. Panagiotis Moulos.  
- The analysis steps presented here are based on his course material and lectures within the MSc program.

---

## References

* HISAT2: https://daehwankimlab.github.io/hisat2/
* SAMtools: http://www.htslib.org/
* metaseqR2: Bioconductor

---
