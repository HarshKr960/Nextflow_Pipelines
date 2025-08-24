**RNA-seq Pipeline (Nextflow DSL2)**

This repository contains a Nextflow DSL2 pipeline for processing RNA-seq paired-end data.
It automates the full workflow from raw FASTQ files to a gene count matrix ready for downstream analysis (e.g., DESeq2, edgeR).

**📌 Workflow Steps**
SetupDirectories → Create output directories (fastqc, trimmed, hisat2, bam, counts).
FastQC → Run quality control on raw FASTQ reads.
Trimmomatic → Adapter removal and quality trimming.
HISAT2 → Alignment of reads to the reference genome.
SamtoolsSort → Convert SAM → BAM and sort alignments.
FeatureCounts → Generate a gene-level count matrix from BAM files.

**Inputs**
  
  Paired-end FASTQ files (*_1.fastq.gz, *_2.fastq.gz)
  
  HISAT2 index (built from reference genome)
  
  GTF annotation file
  
  Adapter file for Trimmomatic

**Outputs**
  
  fastqc/ → Quality reports (.html, .zip)
  
  trimmed/ → Adapter/quality-trimmed reads
  
  hisat2/ → SAM alignment files + HISAT2 mapping stats
 
  bam/ → Sorted BAM files
  
  counts/ → featureCounts.txt (gene counts) + summary
