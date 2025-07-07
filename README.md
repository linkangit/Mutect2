# Somatic Variant Calling with GATK4 Mutect2

This guide walks you through a reproducible workflow for calling somatic variants from whole genome sequencing (WGS) data using [GATK4 Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), following GATK Best Practices.  
All paths are relative to a root directory named `Project`.

---

## Table of Contents

1. [Download Sequencing Reads](#1-download-sequencing-reads)
2. [Subset Reads (Optional)](#2-subset-reads-optional)
3. [Prepare Reference Files](#3-prepare-reference-files)
4. [Quality Control (FastQC)](#4-quality-control-fastqc)
5. [Trim Adapters and Merge Reads](#5-trim-adapters-and-merge-overlapping-reads)
6. [Map Reads to Reference (BWA-MEM)](#6-map-reads-to-reference-bwa-mem)
7. [Mark Duplicates and Sort](#7-mark-duplicates-and-sort)
8. [Base Quality Score Recalibration (BQSR)](#8-base-quality-score-recalibration-bqsr)
9. [Collect Alignment and Insert Size Metrics](#9-collect-alignment-and-insert-size-metrics)
10. [Download Mutect2 Supporting Files](#10-download-mutect2-supporting-files)
11. [Call Somatic Variants (Mutect2)](#11-call-somatic-variants-mutect2)
12. [Estimate Cross-Sample Contamination](#12-estimate-cross-sample-contamination)
13. [Estimate Read Orientation Artifacts](#13-estimate-read-orientation-artifacts)
14. [Filter Mutect2 Calls](#14-filter-mutect2-calls)
15. [Annotate Variants (Funcotator)](#15-annotate-variants-funcotator)
16. [Workflow Summary Table](#workflow-summary-table)

---

## 1. Download Sequencing Reads

First, create a directory for your reads and download the paired-end FASTQ files for both normal and tumor samples.

```bash
mkdir -p Project/reads
cd Project/reads

# Normal (Duodenal tissue)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz

# Tumor (PDAC cell line)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz
```

---

## 2. Subset Reads (Optional)

For demonstration or testing, you can subset a smaller number of reads using `seqtk`.

```bash
seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz 1000000 > HG008-N-D_subset_R1.fastq
seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz 1000000 > HG008-N-D_subset_R2.fastq
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz 1000000 > HG008-T_subset_R1.fastq
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz 1000000 > HG008-T_subset_R2.fastq
```

---

## 3. Prepare Reference Files

Download and prepare the reference genome and known variant sites.

```bash
mkdir -p Project/supporting_files/hg38
cd Project/supporting_files/hg38

# Reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Index and dictionary
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

# Known sites for BQSR
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```

---

## 4. Quality Control (FastQC)

Run FastQC to assess the quality of your reads.

```bash
mkdir -p Project/reads/fastqc

fastqc Project/reads/HG008-N-D_subset_R1.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-N-D_subset_R2.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-T_subset_R1.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-T_subset_R2.fastq -o Project/reads/fastqc
```

---

## 5. Trim Adapters and Merge Overlapping Reads

Trim adapters with Trimmomatic and merge overlapping reads with BBMerge.

```bash
mkdir -p Project/reads/trimmed

# Download adapter file if needed
wget -O Project/reads/TruSeq3-PE.fa https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

# Adapter trimming (Normal)
trimmomatic PE -threads 4 \
  Project/reads/HG008-N-D_subset_R1.fastq Project/reads/HG008-N-D_subset_R2.fastq \
  Project/reads/trimmed/HG008-N-D_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-N-D_trimmed_R1_unpaired.fastq \
  Project/reads/trimmed/HG008-N-D_trimmed_R2_paired.fastq Project/reads/trimmed/HG008-N-D_trimmed_R2_unpaired.fastq \
  ILLUMINACLIP:Project/reads/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Adapter trimming (Tumor)
trimmomatic PE -threads 4 \
  Project/reads/HG008-T_subset_R1.fastq Project/reads/HG008-T_subset_R2.fastq \
  Project/reads/trimmed/HG008-T_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-T_trimmed_R1_unpaired.fastq \
  Project/reads/trimmed/HG008-T_trimmed_R2_paired.fastq Project/reads/trimmed/HG008-T_trimmed_R2_unpaired.fastq \
  ILLUMINACLIP:Project/reads/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Merge overlapping reads
mkdir -p Project/reads/merged

bbmerge.sh in1=Project/reads/trimmed/HG008-N-D_trimmed_R1_paired.fastq \
           in2=Project/reads/trimmed/HG008-N-D_trimmed_R2_paired.fastq \
           out=Project/reads/merged/HG008-N-D_merged.fastq

bbmerge.sh in1=Project/reads/trimmed/HG008-T_trimmed_R1_paired.fastq \
           in2=Project/reads/trimmed/HG008-T_trimmed_R2_paired.fastq \
           out=Project/reads/merged/HG008-T_merged.fastq
```

---

## 6. Map Reads to Reference (BWA-MEM)

Align the trimmed reads to the reference genome using BWA-MEM.

```bash
mkdir -p Project/aligned

# Index reference for BWA
bwa index Project/supporting_files/hg38/hg38.fa

# Align reads
bwa mem -t 4 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" Project/supporting_files/hg38/hg38.fa Project/reads/trimmed/HG008-N-D_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-N-D_trimmed_R2_paired.fastq > Project/aligned/HG008-N-D.paired.sam

bwa mem -t 4 -R "@RG\tID:HG008-T\tPL:ILLUMINA\tSM:HG008-T" Project/supporting_files/hg38/hg38.fa Project/reads/trimmed/HG008-T_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-T_trimmed_R2_paired.fastq > Project/aligned/HG008-T.paired.sam
```

---

## 7. Mark Duplicates and Sort

Use GATK to mark duplicate reads and sort the alignments.

```bash
gatk MarkDuplicatesSpark -I Project/aligned/HG008-N-D.paired.sam -O Project/aligned/HG008-N-D_sorted_dedup_reads.bam
gatk MarkDuplicatesSpark -I Project/aligned/HG008-T.paired.sam -O Project/aligned/HG008-T_sorted_dedup_reads.bam
```

---

## 8. Base Quality Score Recalibration (BQSR)

Perform base quality score recalibration using known variant sites.

```bash
# Build recalibration model
gatk BaseRecalibrator -I Project/aligned/HG008-N-D_sorted_dedup_reads.bam -R Project/supporting_files/hg38/hg38.fa --known-sites Project/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O Project/aligned/HG008-N-D_recal_data.table

gatk BaseRecalibrator -I Project/aligned/HG008-T_sorted_dedup_reads.bam -R Project/supporting_files/hg38/hg38.fa --known-sites Project/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf -O Project/aligned/HG008-T_recal_data.table

# Apply recalibration
gatk ApplyBQSR -I Project/aligned/HG008-N-D_sorted_dedup_reads.bam -R Project/supporting_files/hg38/hg38.fa --bqsr-recal-file Project/aligned/HG008-N-D_recal_data.table -O Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam

gatk ApplyBQSR -I Project/aligned/HG008-T_sorted_dedup_reads.bam -R Project/supporting_files/hg38/hg38.fa --bqsr-recal-file Project/aligned/HG008-T_recal_data.table -O Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam
```

---

## 9. Collect Alignment and Insert Size Metrics

Collect metrics to assess alignment quality and insert size distribution.

```bash
gatk CollectAlignmentSummaryMetrics R=Project/supporting_files/hg38/hg38.fa I=Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam O=Project/aligned/HG008-N-D_alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam OUTPUT=Project/aligned/HG008-N-D_insert_size_metrics.txt HISTOGRAM_FILE=Project/aligned/HG008-N-D_insert_size_histogram.pdf

gatk CollectAlignmentSummaryMetrics R=Project/supporting_files/hg38/hg38.fa I=Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam O=Project/aligned/HG008-T_alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam OUTPUT=Project/aligned/HG008-T_insert_size_metrics.txt HISTOGRAM_FILE=Project/aligned/HG008-T_insert_size_histogram.pdf
```

---

## 10. Download Mutect2 Supporting Files

Download additional resources required for somatic variant calling.

```bash
mkdir -p Project/supporting_files/mutect2_supporting_files
cd Project/supporting_files/mutect2_supporting_files

# gnomAD resource
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

# Panel of Normals (PoN)
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

# Exome calling intervals
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list
```

---

## 11. Call Somatic Variants (Mutect2)

Run Mutect2 to call somatic variants.

```bash
mkdir -p Project/results

gatk Mutect2 \
  -R Project/supporting_files/hg38/hg38.fa \
  -I Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  -I Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  -tumor HG008-T \
  -normal HG008-N-D \
  --germline-resource Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals Project/supporting_files/mutect2_supporting_files/1000g_pon.hg38.vcf.gz \
  -O Project/results/HG008_somatic_variants_mutect2.vcf.gz \
  --f1r2-tar-gz Project/results/HG008_f1r2.tar.gz
```

---

## 12. Estimate Cross-Sample Contamination

Estimate contamination using `GetPileupSummaries` and `CalculateContamination`.

```bash
# Tumor sample
gatk GetPileupSummaries \
  -I Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  -V Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  -L Project/supporting_files/mutect2_supporting_files/exome_calling_regions.v1.1.interval_list \
  -O Project/results/HG008-T_getpileupsummaries.table

# Normal sample
gatk GetPileupSummaries \
  -I Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  -V Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  -L Project/supporting_files/mutect2_supporting_files/exome_calling_regions.v1.1.interval_list \
  -O Project/results/HG008-N-D_getpileupsummaries.table

# Calculate contamination
gatk CalculateContamination \
  -I Project/results/HG008-T_getpileupsummaries.table \
  -matched Project/results/HG008-N-D_getpileupsummaries.table \
  -O Project/results/HG008_pair_calculatecontamination.table
```

---

## 13. Estimate Read Orientation Artifacts

Learn the read orientation model for filtering artifacts.

```bash
gatk LearnReadOrientationModel \
  -I Project/results/HG008_f1r2.tar.gz \
  -O Project/results/read-orientation-model.tar.gz
```

---

## 14. Filter Mutect2 Calls

Filter the raw somatic variant calls.

```bash
gatk FilterMutectCalls \
  -V Project/results/HG008_somatic_variants_mutect2.vcf.gz \
  -R Project/supporting_files/hg38/hg38.fa \
  --contamination-table Project/results/HG008_pair_calculatecontamination.table \
  --ob-priors Project/results/read-orientation-model.tar.gz \
  -O Project/results/HG008_somatic_variants_filtered_mutect2.vcf
```

---

## 15. Annotate Variants (Funcotator)

Annotate filtered variants using GATK Funcotator.

```bash
gatk Funcotator \
  --variant Project/results/HG008_somatic_variants_filtered_mutect2.vcf \
  --reference Project/supporting_files/hg38/hg38.fa \
  --ref-version hg38 \
  --data-sources-path Project/tools/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.8.hg38.20230908s \
  --output Project/results/HG008_somatic_variants_functotated.vcf \
  --output-file-format VCF
```

---

## Workflow Summary Table

| Step | Purpose                                   | Output Directory/Example                                 |
|------|-------------------------------------------|---------------------------------------------------------|
| 1    | Download reads                            | `Project/reads/`                                        |
| 2    | Subset reads (optional)                   | `Project/reads/`                                        |
| 3    | Prepare reference files                   | `Project/supporting_files/hg38/`                        |
| 4    | Quality control (FastQC)                  | `Project/reads/fastqc/`                                 |
| 5    | Trim adapters & merge overlapping reads   | `Project/reads/trimmed/`, `Project/reads/merged/`       |
| 6    | Map reads (BWA-MEM)                       | `Project/aligned/`                                      |
| 7    | Mark duplicates & sort                    | `Project/aligned/`                                      |
| 8    | Base quality recalibration                | `Project/aligned/`                                      |
| 9    | Collect alignment/insert size metrics     | `Project/aligned/`                                      |
| 10   | Download Mutect2 supporting files         | `Project/supporting_files/mutect2_supporting_files/`     |
| 11   | Call somatic variants (Mutect2)           | `Project/results/`                                      |
| 12   | Estimate cross-sample contamination       | `Project/results/`                                      |
| 13   | Estimate read orientation artifacts       | `Project/results/`                                      |
| 14   | Filter Mutect2 calls                      | `Project/results/`                                      |
| 15   | Annotate variants (Funcotator)            | `Project/results/`                                      |

---

## Notes

- This workflow assumes all required tools (GATK4, BWA, Samtools, FastQC, Trimmomatic, BBMerge, seqtk) are installed and available in your `PATH`.
- Adjust thread counts and file paths as needed for your environment.
- For more details on each tool, refer to their official documentation.
